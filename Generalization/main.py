import os
import sys
import sympy as sp
import yaml_loader as yl
import validation as valid
from replace import fix_exp
from sympy.printing.cxx import CXX11CodePrinter
from mako.template import Template
from mako.lookup import TemplateLookup
import emit_header as em
from sympy.printing.cxx import cxxcode
# for mako look at templates
lookup = TemplateLookup(directories=["templates"])


class TT(sp.Function):
    """
    TT(n, x) = sum_{k=0..n} x^k / k!
    Türev: d/dx TT(n,x) = TT(n-1,x)  (n>=1), TT(0,x)' = 0
    """
    @classmethod
    def eval(cls, n, x):
        # İstersen n küçükse burada expand edebilirsin; şimdilik dokunmuyoruz.
        return None

    def fdiff(self, argindex=2):
        # argindex=2 -> x'e göre türev
        if argindex != 2:
            raise ValueError("TT sadece x'e göre türev destekliyor.")
        n, x = self.args
        if n.is_Integer:
            if int(n) <= 0:
                return sp.Integer(0)
            return TT(n - 1, x)
        # n sembolikse: konservatif kal
        return sp.Function.fdiff(self, argindex)


def emit_expr(expr, add_dispersion):
    if add_dispersion:
        # Krypton: pow istemiyoruz
        return cxxcode(expr)
    else:
        # LJ, Mie, vb.
        return ccode_fastpow(expr)

class TempAllocator:
    def __init__(self):
        self.idx = 0

    def new(self):
        name = f"x{self.idx}"
        self.idx += 1
        return name

tmp_alloc = TempAllocator()



class FastPowPrinter(CXX11CodePrinter):
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        if exp.is_Integer and -20 <= int(exp) <= 20:
            return f"fast_pow({self._print(base)}, {int(exp)})"
        return f"std::pow({self._print(base)}, {self._print(exp)})"


_fastpow_printer = FastPowPrinter()

def ccode_fastpow(e):
    return _fastpow_printer.doprint(e)

def lower_tt_series(expr):
    """
    Lowers truncated exponential series TT(n, x) into
    recurrence-based numerical computation.
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
        f"double tt_arr[{nmax+1}];",
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

def calculate_force(expr_str, param_names, add_dispersion):
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)

    sym_locals = {
        "r": r,
        "inv_r": inv_r,
    }
     
    sym_locals["TT"] = TT
    
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
        F = eliminate_r(F)

    F = F.doit()   
    # lowering TT series
    F, tt_prelude = lower_tt_series(F)
    
    repl, exprs = sp.cse(F, optimizations="basic")
    reduced = exprs[0]


    temp_lines = [
        f"const double {s} = {fix_exp(emit_expr(rhs, add_dispersion))};"
        for s, rhs in repl
    ]

    temps_code = tt_prelude
    if tt_prelude and temp_lines:
        temps_code += "\n        "
    temps_code += "\n        ".join(temp_lines)

    force_expr = fix_exp(emit_expr(reduced, add_dispersion))

    return temps_code, force_expr

def eliminate_r(expr):
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)
    return sp.simplify(expr.subs(r, 1 / inv_r))


def generate_soa(cfg: dict, out_dir: str):
    #  for soa.h per potential the parameters and potential name
    classname = cfg["classname"]
    params = cfg["parameters"]
    param_names = list(params.keys())
    # handling cutoff
    cutoff = cfg.get("cutoff")
    cutoff_enabled = cutoff is not None
    cutoff_value = float(cutoff) if cutoff_enabled else 0.0
    # check for dispersion for Krypton
    add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])
    # temporary code extraction and force expression
    temps_code, force_expr = calculate_force(
        cfg["expr_str"], param_names, add_dispersion
    )
    #
    header_tpl = lookup.get_template("soa_functor.h.mako")
    uses_mass = ("p1m" in cfg["expr_str"]) or ("p2m" in cfg["expr_str"])

    ctx = {
        "classname": classname,
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

    hpath = os.path.join(out_dir, f"{classname}_SoA.h")
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

    for cfg in cfgs:
        classname = cfg["classname"]
        params = cfg["parameters"]
        param_names = list(params.keys())

        add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])

        temps, force = calculate_force(
            cfg["expr_str"],
            param_names,
            add_dispersion,
        )

        cpp = em.emit_header(
            classname,
            temps,
            force,
            cfg["newton3_default"],
            cfg["eps_guard"],
            param_names,
            add_dispersion,
        )

        out_path = os.path.join(out_dir, cfg["filename"])
        with open(out_path, "w") as f:
            f.write(cpp)

        if cfg.get("generate_soa", False):
            generate_soa(cfg, out_dir)

        print("[LEVEL-2] Generated:", cfg["filename"])

if __name__ == "__main__":
    main()
