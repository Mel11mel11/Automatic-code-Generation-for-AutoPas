
import re
from typing import Iterable

# in need can be expanded
FUNCS: list[str] = [
    "pow", "sqrt", "cbrt", "fabs",
    "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "sinh", "cosh", "tanh",
    "exp", "log", "log10",
    "floor", "ceil", "fmod", "hypot",
    "erf", "erfc"
]

def make_pattern(func: str) -> re.Pattern:
    return re.compile(rf"(?<!std::)\b{func}\s*\(", re.UNICODE)
_POW_INT = re.compile(
    r'\b(?:std::\s*)?pow\s*\(\s*([^,()]+|\([^()]*\))\s*,\s*(-?\d+)\s*\)'
)

def replace_pow(expression: str) -> str:
    def repl(m):
        base = m.group(1).strip()
        exp  = int(m.group(2))
        if -12 <= exp <= 13:
            if exp == 0:  return "1.0"
            if exp == 1:  return base
            return f"fast_pow({base}, {exp})"
        else:
            return f"std::pow({base}, {exp})"
    out = _POW_INT.sub(repl, expression)
    out = re.sub(r'\bpow\s*\(', 'std::pow(', out)
    return out

def fix_exp(expr: str, funcs: Iterable[str] = FUNCS) -> str:
    out = replace_pow(expr)
    for fn in funcs:
        if fn == "pow":
            continue
        out = make_pattern(fn).sub(f"std::{fn}(", out)
    
    out = re.sub(r'\bstd::\s*std::', 'std::', out)
    return out

