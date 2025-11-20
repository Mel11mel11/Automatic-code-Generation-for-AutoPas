#!/usr/bin/env python3
import os
import sys
import yaml
import sympy as sp
import textwrap
from replace import fix_exp

def _validate_one(data: dict) -> dict:
    name = data.get("name")
    expr = data.get("expression")
    out = data.get("output", {}) or {}
    classname = out.get("classname")
    filename = out.get("filename")
    options = data.get("options", {}) or {}
    params = data.get("parameters", {}) or {}

    if not all([name, expr, classname, filename]):
        raise ValueError("Missing required fields")

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
    with open(path, "r") as f:
        docs = list(yaml.safe_load_all(f))

    items = []
    for doc in docs:
        if not doc:
            continue
        if isinstance(doc, dict) and "potentials" in doc:
            for it in doc["potentials"]:
                items.append(_validate_one(it))
        elif isinstance(doc, dict):
            items.append(_validate_one(doc))

    return items


# ==========================
#  DISPERSION COEFFICIENTS
# ==========================
def compute_dispersion_coeffs_sympy():
    C6 = sp.Symbol("C6", real=True)
    C8 = sp.Symbol("C8", real=True)
    C10 = sp.Symbol("C10", real=True)

    C12 = C6 * (C10 / C8)**3
    C14 = C8 * (C12 / C10)**3
    C16 = C10 * (C14 / C12)**3
    return C12, C14, C16


# ==========================
#  CALCULATE FORCE
# ==========================
def calculate_force(expr_str: str, param_names: list[str], add_dispersion=False):
    r = sp.Symbol("r", positive=True)
    sym_locals = {"r": r}

    # Parametre sembolleri
    for pname in param_names:
        sym_locals[pname] = sp.Symbol(pname, real=True)

    # Krypton / Tang-Toennies ise:
    if add_dispersion:
        C12, C14, C16 = compute_dispersion_coeffs_sympy()
        sym_locals["C12"] = C12
        sym_locals["C14"] = C14
        sym_locals["C16"] = C16

    # Expression → SymPy
    U = sp.sympify(expr_str, locals=sym_locals)

    # Sum(...) kapat
    U = U.doit()
    U = sp.simplify(U)

    # Kuvvet
    F = -sp.diff(U, r)
    F = F.doit()
    F = sp.simplify(F)

    # C++ kodu
    F_c = sp.ccode(F)
    F_std = fix_exp(F_c)

    return F_std

def emit_header(
    classname, F_code, newton3_default, eps_guard, param_names, add_dispersion
):
    # -----------------------------
    # 1) Constructor param listesi
    # -----------------------------
    if param_names:
        ctor_param_sig = ", ".join([f"double {p}" for p in param_names]) + ", "
    else:
        ctor_param_sig = ""

    ctor_param_sig += f"bool newton3 = {'true' if newton3_default else 'false'}"

    # -----------------------------
    # 2) init-list (newton3 EN SONDA)
    # -----------------------------
    init_list_elems = []

    # önce tüm parametre üyeleri:
    for p in param_names:
        init_list_elems.append(f"_{p}({p})")

    # Krypton için dispersion:
    if add_dispersion:
        # ctor body içinde dolduruyoruz
        init_list_elems += ["_C12(0)", "_C14(0)", "_C16(0)"]

    # en sonda newton3:
    init_list_elems.append("_newton3(newton3)")

    init_list = ", ".join(init_list_elems)

    # -----------------------------
    # 3) Parameter aliases (AoSFunctor içinde)
    # -----------------------------
    alias_list = [f"const double {p} = _{p};" for p in param_names]

    # Mie potansiyeli için C alias
    if classname.lower().startswith("mie"):
        alias_list.append("const double C = _C;")

    # Krypton için aliaslar
    if add_dispersion:
        alias_list += [
            "const double C12 = _C12;",
            "const double C14 = _C14;",
            "const double C16 = _C16;",
        ]

    local_aliases = "\n        ".join(alias_list)

    # -----------------------------
    # 4) member declarations (newton3 EN SONDA)
    # -----------------------------
    member_decls_list = []

    # önce tüm parametreler
    for p in param_names:
        member_decls_list.append(f"double _{p};")

    # Mie potansiyeli için C değişkeni
    if classname.lower().startswith("mie"):
        member_decls_list.append("double _C;")

    # Krypton için yeni C değerleri
    if add_dispersion:
        member_decls_list += ["double _C12;", "double _C14;", "double _C16;"]

    # newton3 en sonunda
    member_decls_list.append("bool _newton3;")

    # İşte sende eksik olan satır:
    member_decls = "\n    ".join(member_decls_list)

    # -----------------------------
    # 5) Mie runtime C hesaplama
    # -----------------------------
    if classname.lower().startswith("mie"):
        mie_runtime = """
        _C = (_n / (_n - _m)) * std::pow(_n / _m, _m / (_n - _m));
        """
    else:
        mie_runtime = ""

    # -----------------------------
    # 6) Krypton runtime hesaplama
    # -----------------------------
    if add_dispersion:
        dispersion_runtime = """
        _C12 = _C6 * std::pow(_C10 / _C8, 3.0);
        _C14 = _C8 * std::pow(_C12 / _C10, 3.0);
        _C16 = _C10 * std::pow(_C14 / _C12, 3.0);
        """
    else:
        dispersion_runtime = ""

    # -----------------------------
    # 7) Final C++ header
    # -----------------------------
    body = f"""
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <class Particle_T>
class {classname} : public Functor<Particle_T> {{
public:
    explicit {classname}({ctor_param_sig})
        : {init_list}
    {{
        {mie_runtime}
        {dispersion_runtime}
    }}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {{
        const auto& ra = p1.getR();
        const auto& rb = p2.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        constexpr double EPS = {eps_guard};
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;

        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        // Parameter aliases
        {local_aliases}
        

        
        const double p1m = p1.getMass();
        const double p2m = p2.getMass();

        const double Fmag = {F_code};

        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double,3> F{{fx, fy, fz}};
        p1.addF(F);
        if (_newton3) {{
            p2.subF(F);
        }}
    }}

    bool allowsNewton3() const {{ return true; }}
    bool usesNewton3() const {{ return _newton3; }}

private:
    {member_decls}
}};
"""
    return textwrap.dedent(body)


# ==========================
#  MAIN
# ==========================
def main():
    if len(sys.argv) < 3:
        print("Usage: python3 main.py <spec.yaml> <out_dir/>")
        sys.exit(1)

    yaml_path = sys.argv[1]
    out_dir = sys.argv[2]

    os.makedirs(out_dir, exist_ok=True)

    cfg_list = load_yaml_many(yaml_path)

    for cfg in cfg_list:
        name = cfg["name"]
        params = cfg["parameters"]
        param_names = list(params.keys())

        # Krypton/Tang–Toennies kontrolü
        add_dispersion = all(x in param_names for x in ["C6", "C8", "C10"])

        F_std = calculate_force(cfg["expr_str"], param_names, add_dispersion)

        header = emit_header(
            classname=cfg["classname"],
            F_code=F_std,
            newton3_default=cfg["newton3_default"],
            eps_guard=cfg["eps_guard"],
            param_names=param_names,
            add_dispersion=add_dispersion
        )

        with open(os.path.join(out_dir, cfg["filename"]), "w") as f:
            f.write(header)

        print(f"[ok] Functor created: {cfg['filename']}")


if __name__ == "__main__":
    main()
