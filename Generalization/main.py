import os
import sys
import yaml_loader as yl 
import sympy as sp
import textwrap
import validation as valid
from replace import fix_exp
from sympy.printing.cxx import CXX11CodePrinter
from mako.template import Template
from mako.lookup import TemplateLookup
 # for generate_soa
import emit_header as em  # for emit_header
lookup = TemplateLookup(directories=["templates"])

def generate_soa(cfg: dict, out_dir: str):

    classname = cfg["classname"]
    params = cfg["parameters"]
    param_names = list(params.keys())

    # Force expression
    add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])

    temps_code, force_expr = calculate_force(cfg["expr_str"], param_names, add_dispersion)

    # Load templates
    header_tpl = lookup.get_template("soa_functor.h.mako")
    uses_mass = ("p1m" in cfg["expr_str"]) or ("p2m" in cfg["expr_str"])


    # Context for templates
    ctx = {
        "classname": classname,
        "parameters": params,
        "temps_code": temps_code,
        "force_expr": force_expr,
        "newton3_default": cfg["newton3_default"],
        "eps_guard": cfg["eps_guard"],
        "uses_mass": uses_mass, 
    }

    # Render
    header_code = header_tpl.render(**ctx)
    

    # Output
    hpath = os.path.join(out_dir, f"{classname}_SoA.h")
   

    with open(hpath, "w") as f:
        f.write(header_code)
    

    print("[SoA] Generated:", hpath)

#class FastPowPrinter(CXX11CodePrinter):
   # def _print_Pow(self, expr):
       # base, exp = expr.as_base_exp()
        #if exp.is_Integer and 0 <= exp <= 20:
           # return f"fast_pow({self._print(base)}, {int(exp)})"
        #return f"std::pow({self._print(base)}, {self._print(exp)})"

def eliminate_r(expr):
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)

    # r  -> 1/inv_r
    expr = expr.subs(r, 1/inv_r)

    # sonra sadeleştir (çok önemli)
    expr = sp.simplify(expr)

    return expr

def eliminate_all_r(expr):
    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)

    def repl(node):
        if isinstance(node, sp.Pow):
            base, exp = node.as_base_exp()
            if base.is_Symbol and base.name == "r" and exp.is_integer:
                n = int(exp)
                if n > 0:
                    return 1 / (inv_r ** n)
                else:
                    return inv_r ** (-n)
        return node

    expr = expr.replace(lambda x: isinstance(x, sp.Pow), repl)
    expr = expr.subs(r, 1/inv_r)
    return sp.simplify(expr)

class FastPowPrinter(CXX11CodePrinter):
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()

        # integer exponent
        if exp.is_Integer:
            e = int(exp)

            # small positive integers -> fast_pow
            if 0 <= e <= 20:
                return f"fast_pow({self._print(base)}, {e})"

            # negative integer exponent -> 1 / fast_pow(base, -e)  (fast_pow only needs positive)
            if -20 <= e < 0:
                return f"(1.0/fast_pow({self._print(base)}, {-e}))"

        # fallback
        return f"std::pow({self._print(base)}, {self._print(exp)})"

    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        if exp.is_Integer and 0 <= int(exp) <= 20:
            return f"fast_pow({self._print(base)}, {int(exp)})"
        return f"std::pow({self._print(base)}, {self._print(exp)})"
_fastpow_printer = FastPowPrinter()

def ccode_fastpow(e):
    return _fastpow_printer.doprint(e)

def optimize_level2(expr):

    expr = sp.simplify(expr)

    r = sp.Symbol("r", positive=True)
    inv_r = sp.Symbol("inv_r", positive=True)

    # ----- 1/r^n → inv_r^n -----
    def conv(node):
        if isinstance(node, sp.Pow):
            base, exp = node.as_base_exp()
            if base == r and exp.is_integer and exp < 0:
                return inv_r ** int(-exp)
        return node

    expr = expr.replace(lambda x: isinstance(x, sp.Pow), conv)
    expr = sp.simplify(expr)

    # ----- CSE -----
    repl, exprs = sp.cse(expr, optimizations="basic")
    reduced = exprs[0]

    # ----- exp(...) factoring -----
    exp_terms = [t for t in reduced.atoms(sp.exp)]
    exp_map = {}
    for i, t in enumerate(exp_terms):
        sym = sp.Symbol(f"exp_tmp_{i}")
        exp_map[t] = sym

    if exp_map:
        reduced = reduced.xreplace(exp_map)

    return repl, reduced, exp_map
def factor_inv_r(expr):
    inv_r = sp.Symbol("inv_r", positive=True)

    expr = expr.replace(inv_r**-2, 1/(inv_r*inv_r))
    expr = expr.replace(inv_r**-6, 1/(inv_r**6))
    return sp.simplify(expr)

def compute_dispersion_coeffs_sympy():
    C6, C8, C10 = sp.symbols("C6 C8 C10")
    C12 = C6 * (C10 / C8) ** 3
    C14 = C8 * (C12 / C10) ** 3
    C16 = C10 * (C14 / C12) ** 3
    return C12, C14, C16


def calculate_force(expr_str, param_names, add_dispersion):
    r = sp.Symbol("r", positive=True)
    sym_locals = {"r": r}

    # parameters
    for p in param_names:
        sym_locals[p] = sp.Symbol(p)

    # p1m / p2m
    sym_locals["p1m"] = sp.Symbol("p1m")
    sym_locals["p2m"] = sp.Symbol("p2m")

    # C12-C16 if required
    if add_dispersion:
        C12, C14, C16 = compute_dispersion_coeffs_sympy()
        sym_locals["C12"] = C12
        sym_locals["C14"] = C14
        sym_locals["C16"] = C16

    # Build U, compute F = -dU/dr
    U = sp.sympify(expr_str, locals=sym_locals)
    F = -sp.diff(U, r)

    repl, reduced, exp_map = optimize_level2(F)

    #repl, reduced, exp_map = optimize_level2(F)
    reduced = eliminate_all_r(reduced)
    # reduced = eliminate_r(reduced)
    repl = [(s, eliminate_all_r(rhs)) for (s, rhs) in repl]
    reduced = eliminate_all_r(reduced)

    # (opsiyonel ama iyi) tekrar CSE uygula
    repl2, exprs2 = sp.cse(reduced, optimizations="basic")
    reduced = exprs2[0]
    repl = list(repl) + list(repl2)
    temp_lines = []
    for symb, rhs in repl:
        temp_lines.append(
            f"const double {symb} = {fix_exp(ccode_fastpow(rhs))};"
        )

    # exp tmp’s
    exp_lines = []
    for t, sym in exp_map.items():
        arg = t.args[0]
        exp_lines.append(
            f"const double {sym} = std::exp({fix_exp(ccode_fastpow(arg))});"
        )

    temps_code = (
        "\n        ".join(temp_lines + exp_lines)
        if (temp_lines or exp_lines)
        else ""
    )

    force_expr = fix_exp(ccode_fastpow(reduced))

    return temps_code, force_expr

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 main.py <spec.yaml> <out_dir>")
        return

    yaml_path = sys.argv[1]
    out_dir = sys.argv[2] # output in the second argument
    os.makedirs(out_dir, exist_ok=True)

    raw_docs = yl.load_yaml(yaml_path)

    cfgs = []
    
    for doc in raw_docs:
            cfgs.append(valid.validate_one(doc))
   

        
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
            add_dispersion
            )

        cpp = em.emit_header(
            classname,
            temps,
            force,
            cfg["newton3_default"],   # ✔ doğru
            cfg["eps_guard"],
            param_names,
            add_dispersion,
        )

        out_path = os.path.join(out_dir, cfg["filename"])
        with open(out_path, "w") as f:
            f.write(cpp)

    
        if cfg["generate_soa"]:
            generate_soa(cfg, out_dir)

        print("[LEVEL-2] Generated:", cfg["filename"])


if __name__ == "__main__":
    main()
