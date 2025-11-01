#!/usr/bin/env python3
import os
import sys
import yaml
import sympy as sp
import textwrap
from replace import fix_exp  # pow -> std::pow, etc.


def _validate_one(data: dict) -> dict:
    if not isinstance(data, dict):
        raise ValueError("YAML item must be a mapping (dict).")
    name = data.get("name")
    expr = data.get("expression")
    out = data.get("output", {}) or {}
    classname = out.get("classname")
    filename = out.get("filename")
    options = data.get("options", {}) or {}
    params = data.get("parameters", {}) or {}

    if not all([name, expr, classname, filename]):
        raise ValueError(
            "YAML item validation error: 'name', 'expression', 'output.classname', and 'output.filename' are required."
        )

    for k, v in params.items():
        if not isinstance(k, str):
            raise ValueError(f"Parameter name must be string, got: {k!r}")
        if not isinstance(v, (int, float, str, type(None))):
            raise ValueError(f"Parameter '{k}' must be number or string, got: {type(v)}")

    newton3_default = bool(options.get("newton3", True))
    eps_guard = float(options.get("avoid_r2_zero", 1e-24))

    return {
        "name": name,
        "expr_str": expr,
        "classname": classname,
        "filename": filename,
        "newton3_default": newton3_default,
        "eps_guard": eps_guard,
        "parameters": params,
    }


def load_yaml_many(path: str) -> list[dict]:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"YAML not found: {path}")

    with open(path, "r", encoding="utf-8") as f:
        docs = list(yaml.safe_load_all(f))

    items: list[dict] = []
    for doc in docs:
        if doc is None:
            continue
        # case A: potentials: [ {...}, {...} ]
        if isinstance(doc, dict) and "potentials" in doc and isinstance(doc["potentials"], list):
            for it in doc["potentials"]:
                items.append(_validate_one(it))
        # case B: single flat item
        elif isinstance(doc, dict):
            items.append(_validate_one(doc))
        # other types are ignored

    if not items:
        raise ValueError("No valid potentials found. Use '---' multi-doc or 'potentials: [ ... ]' list.")
    return items


def calculate_force(expr_str: str, param_names: list[str]):
    r = sp.Symbol("r", positive=True)
    sym_locals = {"r": r}
    for pname in param_names:
        sym_locals[pname] = sp.Symbol(pname, real=True)

    U = sp.sympify(expr_str, locals=sym_locals)
    F = -sp.diff(U, r)
    F = sp.nsimplify(F)
    F = sp.N(F)

    F_c = sp.ccode(F)
    F_std = fix_exp(F_c)
    return F, F_c, F_std


def emit_header(
    classname: str,
    F_code: str,
    newton3_default: bool,
    eps_guard: float,
    src_yaml: str,
    param_names: list[str],
) -> str:
   

    member_inits = ", ".join([f"_{p}({p})" for p in param_names])
    member_inits = (", " + member_inits) if member_inits else ""
    local_aliases = "\n        ".join([f"const double {p} = _{p};" for p in param_names]) if param_names else "// (no parameters)"
    member_decls = "\n    ".join([f"double _{p};" for p in param_names]) if param_names else ""
    ctor_param_sig = (", " + ", ".join([f"double {p}" for p in param_names])) if param_names else ""

    body = f"""
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>

template <class Particle_T>
class {classname} : public Functor<Particle_T> {{
public:
    explicit {classname}(bool newton3 = {"true" if newton3_default else "false"}{ctor_param_sig})
        : _newton3(newton3){member_inits} {{}}

    bool allowsNewton3() const override {{ return true; }}
    bool usesNewton3()   const override {{ return _newton3; }}

    void AoSFunctor(Particle_T& a, Particle_T& b) override {{
        // Displacement a -> b (keep the same direction convention as reference)
        const auto& ra = a.getR();
        const auto& rb = b.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        // r^2 and r; guard against r -> 0
        constexpr double EPS = {eps_guard};
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;
        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        // Parameter aliases
        {local_aliases}

        // --- codegen: force magnitude F(r, params) ---
        // Fmag = {F_code}
        const double Fmag = {F_code};

        // Vector force: F = Fmag * r̂
        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        a.addF(fx, fy, fz);
        if (_newton3) {{
            b.addF(-fx, -fy, -fz);
        }}
    }}

private:
    bool _newton3;
    {member_decls}
}};
"""
    return textwrap.dedent(body)


def main():
    # Usage: python3 script.py <spec.yaml> <out_dir>
    if len(sys.argv) < 3:
        print("Usage: python3 <script.py> <spec.yaml> <out_dir>")
        sys.exit(1)

    yaml_path = sys.argv[1]
    out_dir = sys.argv[2]

    try:
        os.makedirs(out_dir, exist_ok=True)
        cfg_list = load_yaml_many(yaml_path)

        # Optional: deduplicate filenames/classnames early
        seen_files, seen_classes = set(), set()
        for cfg in cfg_list:
            if cfg["filename"] in seen_files:
                raise ValueError(f"Duplicate output filename: {cfg['filename']}")
            if cfg["classname"] in seen_classes:
                raise ValueError(f"Duplicate classname: {cfg['classname']}")
            seen_files.add(cfg["filename"])
            seen_classes.add(cfg["classname"])

        for cfg in cfg_list:
            param_names = list(cfg["parameters"].keys())
            _, _, F_std = calculate_force(cfg["expr_str"], param_names)

            header_str = emit_header(
                classname=cfg["classname"],
                F_code=F_std,
                newton3_default=cfg["newton3_default"],
                eps_guard=cfg["eps_guard"],
                src_yaml=yaml_path,
                param_names=param_names,
            )

            out_path = os.path.join(out_dir, cfg["filename"])
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(header_str)

            print(f"[ok] Functor created: {out_path}")

    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(2)

if __name__ == "__main__":
    main()
