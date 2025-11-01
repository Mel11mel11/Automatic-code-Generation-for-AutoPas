# codegen/cxxify.py
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

def fix_exp(expr: str, funcs: Iterable[str] = FUNCS) -> str:
    # makes C++ code std::namespace 
    out = expr
    for fn in funcs:
        out = make_pattern(fn).sub(f"std::{fn}(", out)
    return out
