import os
import sys
import sympy as sp
from dataclasses import dataclass

import yaml_loader as yl  # my helper: reads YAML specs (potentials, params, options)
import validation as valid  # my helper: checks the YAML fields and fills defaults
from replace import fix_exp  # small helper: converts SymPy exp(...) to C++ std::exp(...) etc.
from replace import replace_pow  # replaces std::pow(a, n) by fast_pow(a, n) for small integer n
from sympy.printing.cxx import CXX11CodePrinter  # SymPy printer base for C++ code generation
from mako.lookup import TemplateLookup  # templating for C++ header generation
import emit_header as em  # my helper: writes the final AoS header (C++ functor)
from sympy.printing.cxx import cxxcode  # standard SymPy -> C++ printer

# We use mako templates for SoA functors (templates folder is in the repo)
lookup = TemplateLookup(directories=["templates"])


@dataclass()  # dataclass is just a compact way to store config values
class Optimizations:
    """
    This class stores which symbolic optimizations should be applied.
    I use it to generate several versions and compare performance.
    """

    # after substituting r -> 1/inv_r, optionally simplify the symbolic expression
    simplify_elim: bool = False

    # SymPy CSE = common subexpression elimination (introduces temporaries)
    cse: bool = False
    # basic|expr2 = sp.powsimp(expr)|None
    # repl, exprs = sp.cse(expr2, optimizations="basic")
    # use fast_pow(a, n) for small integer n instead of std::pow(a, n)
    fast_pow: bool = False

    # Note: aggressive_factor is currently disabled in the class on purpose
    # (I tested it, but it is experimental and not stable for all cases)
    # inv_r^3 * inv_r^5  →  inv_r^8
    # aggressive_factor: bool = False


def opt_suffix(opt: Optimizations) -> str:
    """
    Create a short suffix so each generated file/class has a unique name,
    depending on which optimizations are enabled.
    """
    return f"Opt{int(opt.simplify_elim)}{int(opt.cse)}{int(opt.fast_pow)}"


# ---- TT: truncated Taylor series for exp(x) ---------------------------------
# WARNING: In the current version this is always enabled (mandatory),
# because calculate_force() always calls lower_tt_series().

class TT(sp.Function):
    """
    TT(n, x) = sum_{k=0..n} x^k / k!
    I use this as a symbolic placeholder so SymPy can differentiate it.

    Derivative w.r.t. x:
      d/dx TT(n, x) = TT(n-1, x)  for n >= 1
      d/dx TT(0, x) = 0
    """

    @classmethod
    def eval(cls, n, x):
        # Returning None tells SymPy: "do not simplify / evaluate this automatically"
        # We want to keep TT(...) symbolic until we explicitly lower it later.
        return None

    def fdiff(self, argindex=2):
        # fdiff defines the derivative rule for custom SymPy functions
        if argindex != 2:
            raise ValueError("TT sadece x'e göre türev destekliyor.")  # only d/dx is supported

        n, x = self.args

        # If n is an integer, we can apply our closed-form derivative rule
        if n.is_Integer:
            if int(n) <= 0:
                return sp.Integer(0)
            return TT(n - 1, x)

        # If SymPy cannot prove n is integer, fall back to generic behavior
        return sp.Function.fdiff(self, argindex)


# ---- FastPow printer --------------------------------------------------------
class StdPowPrinter(CXX11CodePrinter):
    """
    Strict printer: always prints Pow(base, exp) as std::pow(base, exp).
    This guarantees that fast_pow never appears when opt.fast_pow is False.
    """
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        return f"std::pow({self._print(base)}, {self._print(exp)})"


_stdpow_printer = StdPowPrinter()

def ccode_stdpow(e):
    return _stdpow_printer.doprint(e)

class FastPowPrinter(CXX11CodePrinter):
    """
    Custom C++ printer: prints Pow(base, exp) as fast_pow(base, exp)
    for small integer exponents, otherwise prints std::pow(base, exp).
    """

    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()

        # For small integer exponents we prefer fast_pow (in our C++ codebase)
        if exp.is_Integer and -20 <= int(exp) <= 20:
            return f"fast_pow({self._print(base)}, {int(exp)})"

        # Otherwise use standard pow (works for real exponents etc.)
        return f"std::pow({self._print(base)}, {self._print(exp)})"


_fastpow_printer = FastPowPrinter()


def ccode_fastpow(e):
    """Convenience wrapper around our custom printer."""
    return _fastpow_printer.doprint(e)



def emit_expr(expr, add_dispersion: bool, opt: Optimizations) -> str:
    if add_dispersion:
        return cxxcode(expr)

    if opt.fast_pow:
        return ccode_fastpow(expr)

    return ccode_stdpow(expr)



def lower_tt_series(expr):
    """
    Replace TT(n, x) nodes inside the expression by numeric variables ttN
    and generate a C++ prelude that computes them using a recurrence.

    Important: This is currently applied unconditionally in calculate_force().
    If the expression has no TT(...) terms, it returns (expr, "").
    """
    # Find all TT(...) occurrences in yaml expression
    tts = list(expr.atoms(TT))
    if not tts:
        return expr, ""

    # We currently support only one shared x across all TT(n, x) in the expression
    xs = {tt.args[1] for tt in tts}
    if len(xs) != 1:
        return expr, ""  # not supported -> do nothing

    x = list(xs)[0]

    # Collect the truncation orders that appear
    ns = sorted({int(tt.args[0]) for tt in tts if tt.args[0].is_Integer})
    if not ns:
        return expr, ""  # if n is not integer -> do nothing

    nmax = max(ns)

    # Replace TT(n, x) by symbols ttN inside the SymPy expression
    repl = {TT(sp.Integer(n), x): sp.Symbol(f"tt{n}") for n in ns}
    expr2 = expr.xreplace(repl)

    # Convert x to C++ string
    x_cpp = sp.ccode(x)

    # Generate C++ code that builds TT values iteratively:
    # tt_arr[k] = sum_{i=0..k} x^i / i!
    prelude = [
        f"const double tt_x = {x_cpp};",
        # NOTE: exp_m_tt_x was removed because it was unused
        f"double tt_arr[{nmax + 1}];",
        "tt_arr[0] = 1.0;",
        "double tt_pow = 1.0;",       # x^k
        "double tt_inv_fact = 1.0;",  # 1/k!

        "double tt_sum = 1.0;",       # running sum
        f"for (int k = 1; k <= {nmax}; ++k) {{",
        "  tt_pow *= tt_x;",
        "  tt_inv_fact /= double(k);",
        "  tt_sum += tt_pow * tt_inv_fact;",
        "  tt_arr[k] = tt_sum;",
        "}",
    ]

    # Expose the specific orders we actually need (tt6, tt8, ...)
    for n in ns:
        prelude.append(f"const double tt{n} = tt_arr[{n}];")

    return expr2, "\n        ".join(prelude)


# ---- Krypton helper: derive higher dispersion coefficients -------------------

def compute_dispersion_coeffs_sympy():
    """
    In the Krypton model, C12/C14/C16 are derived from (C6, C8, C10).
    This function defines them symbolically so they can be used in SymPy.
    """
    C6, C8, C10 = sp.symbols("C6 C8 C10")

    # These are model-specific formulas (based on ratios)
    C12 = C6 * (C10 / C8) ** 3
    C14 = C8 * (C12 / C10) ** 3
    C16 = C10 * (C14 / C12) ** 3

    return C12, C14, C16


# Experimental aggressive factoring was removed/commented on purpose.
# I keep it here as future work, but I do not include it in the experiments.

# def aggressive_factor_before_cse(expr):
#     """
#     Experimental: try to factor/powsimp before CSE to expose more common terms.
#     This can also blow up expressions, so it's not enabled by default.
#     """
#     expr = sp.powsimp(expr, force=True)
#     expr = sp.factor_terms(expr)
def detect_use_mass(cfg: dict) -> bool:
    # Gravity model: masses are mandatory (p1m, p2m)
    # LJ/Mie/Krypton: no masses
    name = (cfg.get("classname") or "").lower()
    return name.startswith("gravity")


# ---- r -> inv_r elimination -------------------------------------------------

def eliminate_r(expr, do_simplify: bool):
    """
    Replace r by 1/inv_r (so we can reuse inv_r powers easily in the generated code).
    Optionally run simplify() after substitution.
    """
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)

    expr2 = expr.subs(r, 1 / inv_r)
    return sp.simplify(expr2) if do_simplify else expr2


# ---- main symbolic force generation ----------------------------------------

def calculate_force(expr_str, param_names, add_dispersion: bool, opt: Optimizations):
    """
    Takes a potential expression U(r) as string and returns:
      - temps_code: C++ lines for temporary variables (CSE + TT prelude)
      - force_expr: final force expression (as C++ string)

    Steps:
      1) parse expression
      2) derive force: F = -dU/dr
      3) apply optional symbolic transformations (r elimination, simplify, CSE)
      4) print expression as C++ code
    """
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)

    # locals for sympify: define which symbols/functions can appear in YAML expression
    sym_locals = {"r": r, "inv_r": inv_r, "TT": TT}

    # add parameters (sigma, epsilon, etc.)
    for p in param_names:
        sym_locals[p] = sp.Symbol(p)

    # mass parameters for gravity-like terms
    sym_locals["p1m"] = sp.Symbol("p1m")
    sym_locals["p2m"] = sp.Symbol("p2m")

    # for Krypton: also provide symbolic C12/C14/C16 derived from C6/C8/C10
    if add_dispersion:
        C12, C14, C16 = compute_dispersion_coeffs_sympy()
        sym_locals["C12"] = C12
        sym_locals["C14"] = C14
        sym_locals["C16"] = C16

    # parse potential energy U from YAML string
    U = sp.sympify(expr_str, locals=sym_locals)

    # derive force magnitude: F = -dU/dr
    F = -sp.diff(U, r)

    # For non-dispersion cases we rewrite r -> inv_r (Krypton keeps r form here)
    if not add_dispersion:
        F = eliminate_r(F, do_simplify=opt.simplify_elim)

    # doit() triggers evaluation of derivatives like TT(...)
    F = F.doit()

    # TT lowering: if there is no TT, it returns a no-op prelude
    F, tt_prelude = lower_tt_series(F)

    # Optional experimental factoring (disabled currently)
    # if opt.aggressive_factor:
    #     F = aggressive_factor_before_cse(F)

    # Optional CSE: introduce temporary variables for repeated subexpressions
    if opt.cse:
        repl, exprs = sp.cse(F, optimizations="basic")
        reduced = exprs[0]
    else:
        repl, reduced = [], F

    # Convert CSE temporaries into C++ lines
    #temp_lines = [
      #  f"const double {s} = {fix_exp(emit_expr(rhs, add_dispersion, opt))};"
       # for s, rhs in repl
    #]
    temp_lines = []
    for s, rhs in repl:
        rhs_cpp = emit_expr(rhs, add_dispersion, opt)
        rhs_cpp = fix_exp(rhs_cpp)
        if opt.fast_pow:
             rhs_cpp = replace_pow(rhs_cpp)
        temp_lines.append(f"const double {s} = {rhs_cpp};")
    # Merge TT prelude + CSE temporaries
    temps_code = tt_prelude
    if tt_prelude and temp_lines:
        temps_code += "\n        "
    temps_code += "\n        ".join(temp_lines)

    # Final force expression in C++ form
    force_expr = emit_expr(reduced, add_dispersion, opt)
    force_expr = fix_exp(force_expr)
    if opt.fast_pow:
        force_expr = replace_pow(force_expr)

    return temps_code, force_expr


# ---- SoA codegen path -------------------------------------------------------

def generate_soa(cfg: dict, out_dir: str, opt: Optimizations):
    """
    Generate a SoA functor header using mako template.
    Cutoff is handled as a configuration value passed into the template.
    (The actual cutoff check is expected to be in the generated C++.)
    """
    classname = cfg["classname"]
    params = cfg["parameters"]
    param_names = list(params.keys())

    # cutoff is optional in YAML
    cutoff = cfg.get("cutoff")
    cutoff_enabled = cutoff is not None
    cutoff_value = float(cutoff) if cutoff_enabled else 0.0

    # Krypton/dispersion detection: if C6,C8,C10 exist, we derive higher terms
    add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])

    # generate (temps, force) strings for C++
    temps_code, force_expr = calculate_force(cfg["expr_str"], param_names, add_dispersion, opt)

    # load mako template for SoA functor
    header_tpl = lookup.get_template("soa_functor.h.mako")

    # detect if expression uses masses (gravity)
    use_mass = detect_use_mass(cfg)

    suffix = opt_suffix(opt)

    # context variables for the template
    ctx = {
        "classname": f"{classname}_{suffix}",
        "parameters": params,
        "temps_code": temps_code,
        "force_expr": force_expr,
        "newton3_default": cfg["newton3_default"],
        "eps_guard": cfg["eps_guard"],
        "use_mass": use_mass,
        "cutoff": cutoff_value,
        "cutoff_enabled": cutoff_enabled,
    }

    # render C++ header from template
    header_code = header_tpl.render(**ctx)

    # write result to output directory
    hpath = os.path.join(out_dir, f"{classname}_{suffix}_SoA.h")
    with open(hpath, "w") as f:
        f.write(header_code)

    print("[SoA] Generated:", hpath)


# ---- CLI entry point --------------------------------------------------------

def main():
    """
    Entry point:
      python3 main.py <spec.yaml> <out_dir>
    Generates AoS + optional SoA headers for all YAML potentials and all optimization variants.
    """
    if len(sys.argv) < 3:
        print("Usage: python3 main.py <spec.yaml> <out_dir>")
        return

    yaml_path = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    # read YAML documents (can be multiple separated by ---)
    raw_docs = yl.load_yaml(yaml_path)
    cfgs = []

    # validate and normalize each potential spec
    for doc in raw_docs:
        if not doc:
            continue
        if isinstance(doc, dict) and "potentials" in doc:
            for p in doc["potentials"]:
                cfgs.append(valid.validate_one(p))
        else:
            cfgs.append(valid.validate_one(doc))

    # list of optimization combinations we want to generate (ablation study)
    # NOTE: aggressive_factor is not part of the dataclass currently,
    # so do not pass a 4th boolean here unless you re-add it.
    opt_list = [
        Optimizations(False, False, False),  # opt000 baseline
        Optimizations(True,  False, False),  # opt100 simplify
        Optimizations(False, True,  False),  # opt010 cse
        Optimizations(False, False, True),   # opt001 fast_pow

        Optimizations(True,  False, True),   # opt101 simplify + fast_pow
        Optimizations(True,  True,  False),  # opt110 simplify + cse
        Optimizations(False, True,  True),   # opt011 cse + fast_pow

        Optimizations(True,  True,  True),   # opt111 full (3 opts)
    ]

    # generate headers for every potential and every optimization setting
    for cfg in cfgs:
        base_classname = cfg["classname"]
        params = cfg["parameters"]
        param_names = list(params.keys())

        add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])

    # >>> BUNU EKLE <<<
        use_mass = use_mass = detect_use_mass(cfg)

        for opt in opt_list:
                suffix = opt_suffix(opt)
                classname = f"{base_classname}_{suffix}"

            # get temps + force C++ strings
                temps, force = calculate_force(cfg["expr_str"], param_names, add_dispersion, opt)

            # AoS header code generation (no template here, uses emit_header helper)
                cpp = em.emit_header(
            classname,
            temps,
            force,
            cfg["newton3_default"],
            cfg["eps_guard"],
            param_names,
            add_dispersion,
            use_fast_pow=opt.fast_pow,
            use_mass=use_mass 
)


            # output filename is based on YAML filename + optimization suffix
                base, ext = os.path.splitext(cfg["filename"])
                out_name = f"{base}_{suffix}{ext}"
                out_path = os.path.join(out_dir, out_name)

            # write AoS header to disk
                with open(out_path, "w") as f:
                    f.write(cpp)

            # optionally also generate SoA header
                if cfg.get("generate_soa", False):
                    generate_soa(cfg, out_dir, opt)

                print(" Generated:", out_name)


if __name__ == "__main__":
    main()
