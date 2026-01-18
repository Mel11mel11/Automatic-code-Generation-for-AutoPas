import os
import sys
import sympy as sp
from dataclasses import dataclass

import yaml_loader as yl # for reading the input
import validation as valid # for validating the input
from replace import fix_exp # to fix exp(...) in generated code
from sympy.printing.cxx import CXX11CodePrinter
from mako.lookup import TemplateLookup
import emit_header as em
from sympy.printing.cxx import cxxcode

# for mako look at templates
lookup = TemplateLookup(directories=["templates"])


# ---------------------------
# Optimization configuration
# ---------------------------
@dataclass(frozen=True)
class OptimizationConfig:
    simplify_elim: bool = False  # sp.simplify after r -> inv_r substitution
    cse: bool = False            # SymPy common subexpression elimination
    fast_pow: bool = False       # fast_pow for small integer exponents (vs std::pow)
    aggressive_factor: bool = False  # EXPERIMENTAL

def opt_suffix(opt: OptimizationConfig) -> str:
    return f"Opt{int(opt.simplify_elim)}{int(opt.cse)}{int(opt.fast_pow)}{int(opt.aggressive_factor)}"



# TT symbolic function
class TT(sp.Function):
    """
    TT(n, x) = sum_{k=0..n} x^k / k!
    Derivative: d/dx TT(n,x) = TT(n-1,x)  (n>=1), TT(0,x)' = 0
    """
    @classmethod
    def eval(cls, n, x):
        return None

    def fdiff(self, argindex=2):
        if argindex != 2:
            raise ValueError("TT sadece x'e göre türev destekliyor.")
        n, x = self.args
        if n.is_Integer:
            if int(n) <= 0:
                return sp.Integer(0)
            return TT(n - 1, x)
        return sp.Function.fdiff(self, argindex)



# FastPow printer 
class FastPowPrinter(CXX11CodePrinter):
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        if exp.is_Integer and -20 <= int(exp) <= 20:
            return f"fast_pow({self._print(base)}, {int(exp)})"
        return f"std::pow({self._print(base)}, {self._print(exp)})"


_fastpow_printer = FastPowPrinter()


def ccode_fastpow(e):
    return _fastpow_printer.doprint(e)


def emit_expr(expr, add_dispersion: bool, opt: OptimizationConfig) -> str:
    """
    - Krypton (add_dispersion=True): use plain CXX code (no fast_pow).
    - Others: choose between fast_pow printer and standard cxxcode printer.
    """
    if add_dispersion:
        return cxxcode(expr)

    if opt.fast_pow:
        return ccode_fastpow(expr)

    # fast_pow disabled: standard printer -> std::pow
    return cxxcode(expr)


# TT lowering (mandatory)

def lower_tt_series(expr):
    """
    Lowers truncated exponential series TT(n, x) into recurrence-based numerical computation.
    Applied unconditionally; for non-Krypton potentials, it is typically a no-op.
    """
    tts = list(expr.atoms(TT))
    if not tts:
        return expr, ""

    xs = {tt.args[1] for tt in tts}
    if len(xs) != 1:
        return expr, ""

    x = list(xs)[0]
    ns = sorted({int(tt.args[0]) for tt in tts if tt.args[0].is_Integer})
    if not ns:
        return expr, ""

    nmax = max(ns)

    repl = {TT(sp.Integer(n), x): sp.Symbol(f"tt{n}") for n in ns}
    expr2 = expr.xreplace(repl)

    x_cpp = sp.ccode(x)

    prelude = [
        f"const double tt_x = {x_cpp};",
        "const double exp_m_tt_x = std::exp(-tt_x);",
        f"double tt_arr[{nmax + 1}];",
        "tt_arr[0] = 1.0;",
        "double tt_pow = 1.0;",
        "double tt_inv_fact = 1.0;",
        "double tt_sum = 1.0;",
        f"for (int k = 1; k <= {nmax}; ++k) {{",
        "  tt_pow *= tt_x;",
        "  tt_inv_fact /= double(k);",
        "  tt_sum += tt_pow * tt_inv_fact;",
        "  tt_arr[k] = tt_sum;",
        "}",
    ]

    for n in ns:
        prelude.append(f"const double tt{n} = tt_arr[{n}];")

    return expr2, "\n        ".join(prelude)


def compute_dispersion_coeffs_sympy():
    C6, C8, C10 = sp.symbols("C6 C8 C10")
    C12 = C6 * (C10 / C8) ** 3
    C14 = C8 * (C12 / C10) ** 3
    C16 = C10 * (C14 / C12) ** 3
    return C12, C14, C16
def aggressive_factor_before_cse(expr):
    """
    Experimental algebraic normalization pass to expose common inverse-distance
    subexpressions before CSE.

    WARNING:
    - May increase expression size
    - May change evaluation order (but not semantics)
    - Intended for experimentation only
    """

    # Try to normalize powers (e.g., inv_r^a * inv_r^b -> inv_r^(a+b))
    expr = sp.powsimp(expr, force=True)

    # Factor out common multiplicative terms
    expr = sp.factor_terms(expr)

    return expr



def eliminate_r(expr, do_simplify: bool):
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)
    expr2 = expr.subs(r, 1 / inv_r)
    return sp.simplify(expr2) if do_simplify else expr2


def calculate_force(expr_str, param_names, add_dispersion: bool, opt: OptimizationConfig):
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)

    sym_locals = {"r": r, "inv_r": inv_r, "TT": TT}

    for p in param_names:
        sym_locals[p] = sp.Symbol(p)

    sym_locals["p1m"] = sp.Symbol("p1m")
    sym_locals["p2m"] = sp.Symbol("p2m")

    if add_dispersion:
        C12, C14, C16 = compute_dispersion_coeffs_sympy()
        sym_locals["C12"] = C12
        sym_locals["C14"] = C14
        sym_locals["C16"] = C16

    U = sp.sympify(expr_str, locals=sym_locals)
    F = -sp.diff(U, r)

    if not add_dispersion:
        F = eliminate_r(F, do_simplify=opt.simplify_elim)

    F = F.doit()

    # Mandatory lowering step (Krypton); no-op otherwise
    F, tt_prelude = lower_tt_series(F)

#   Experimental aggressive factoring
    if opt.aggressive_factor:
     F = aggressive_factor_before_cse(F)

# Optional CSE

    # Optional CSE
    if opt.cse:
        repl, exprs = sp.cse(F, optimizations="basic")
        reduced = exprs[0]
    else:
        repl, reduced = [], F

    temp_lines = [
        f"const double {s} = {fix_exp(emit_expr(rhs, add_dispersion, opt))};"
        for s, rhs in repl
    ]

    temps_code = tt_prelude
    if tt_prelude and temp_lines:
        temps_code += "\n        "
    temps_code += "\n        ".join(temp_lines)

    force_expr = fix_exp(emit_expr(reduced, add_dispersion, opt))
    return temps_code, force_expr


def generate_soa(cfg: dict, out_dir: str, opt: OptimizationConfig):
    classname = cfg["classname"]
    params = cfg["parameters"]
    param_names = list(params.keys())

    cutoff = cfg.get("cutoff")
    cutoff_enabled = cutoff is not None
    cutoff_value = float(cutoff) if cutoff_enabled else 0.0

    add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])

    temps_code, force_expr = calculate_force(cfg["expr_str"], param_names, add_dispersion, opt)

    header_tpl = lookup.get_template("soa_functor.h.mako")
    uses_mass = ("p1m" in cfg["expr_str"]) or ("p2m" in cfg["expr_str"])

    suffix = opt_suffix(opt)

    ctx = {
        "classname": f"{classname}_{suffix}",
        "parameters": params,
        "temps_code": temps_code,
        "force_expr": force_expr,
        "newton3_default": cfg["newton3_default"],
        "eps_guard": cfg["eps_guard"],
        "uses_mass": uses_mass,
        "cutoff": cutoff_value,
        "cutoff_enabled": cutoff_enabled,
    }

    header_code = header_tpl.render(**ctx)

    hpath = os.path.join(out_dir, f"{classname}_{suffix}_SoA.h")
    with open(hpath, "w") as f:
        f.write(header_code)

    print("[SoA] Generated:", hpath)


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 main.py <spec.yaml> <out_dir>")
        return

    yaml_path = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    raw_docs = yl.load_yaml(yaml_path)
    cfgs = []

    for doc in raw_docs:
        if not doc:
            continue
        if isinstance(doc, dict) and "potentials" in doc:
            for p in doc["potentials"]:
                cfgs.append(valid.validate_one(p))
        else:
            cfgs.append(valid.validate_one(doc))

    # Minimal set for ablation (5 variants)
    opt_list = [
    OptimizationConfig(False, False, False, False),  # O0000 baseline

    OptimizationConfig(True,  False, False, False),  # O1000 simplify
    OptimizationConfig(False, True,  False, False),  # O0100 CSE
    OptimizationConfig(False, False, True,  False),  # O0010 fast_pow

    OptimizationConfig(True,  False, True,  False),  # O1010 simplify + fast_pow
    OptimizationConfig(True,  True,  False, False),  # O1100 simplify + CSE
    OptimizationConfig(False, True,  True,  False),  # O0110 CSE + fast_pow

    OptimizationConfig(True,  True,  True,  False),  # O1110 full (3 optimizations)

    # --- Experimental variant ---
    OptimizationConfig(False, True,  False, True),   # O0101 CSE + aggressive_factor (EXPERIMENTAL)
]

   

    for cfg in cfgs:
        base_classname = cfg["classname"]
        params = cfg["parameters"]
        param_names = list(params.keys())

        add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])

        for opt in opt_list:
            suffix = opt_suffix(opt)
            classname = f"{base_classname}_{suffix}"

            temps, force = calculate_force(cfg["expr_str"], param_names, add_dispersion, opt)

            cpp = em.emit_header(
                classname,
                temps,
                force,
                cfg["newton3_default"],
                cfg["eps_guard"],
                param_names,
                add_dispersion,
            )

            base, ext = os.path.splitext(cfg["filename"])
            out_name = f"{base}_{suffix}{ext}"
            out_path = os.path.join(out_dir, out_name)

            with open(out_path, "w") as f:
                f.write(cpp)

            if cfg.get("generate_soa", False):
                generate_soa(cfg, out_dir, opt)

            print("[LEVEL-2] Generated:", out_name)


if __name__ == "__main__":
    main()
