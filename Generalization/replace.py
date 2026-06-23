import re
from typing import Iterable

# Math functions that SymPy may print without the std:: namespace prefix
FUNCS: list[str] = [
    "pow", "sqrt", "cbrt", "fabs",
    "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "sinh", "cosh", "tanh",
    "exp", "log", "log10",
    "floor", "ceil", "fmod", "hypot",
    "erf", "erfc",
]


def make_pattern(func: str) -> re.Pattern:
    # Matches the function name only when it is not already preceded by "std::"
    return re.compile(rf"(?<!std::)\b{func}\s*\(", re.UNICODE)


# Matches std::pow or bare pow with an integer exponent so we can swap in fast_pow
_POW_INT = re.compile(
    r'\b(?:std::\s*)?pow\s*\(\s*([^,()]+|\([^()]*\))\s*,\s*(-?\d+)\s*\)'
)


def replace_pow(expression: str) -> str:
    """
    Replace std::pow(base, n) with fast_pow(base, n) for small non-negative
    integer exponents (0–20). Negative or large exponents keep std::pow.
    """
    def repl(m):
        base = m.group(1).strip()
        exp = int(m.group(2))
        if exp == 0:
            return "1.0"
        if exp == 1:
            return base
        if 2 <= exp <= 20:
            return f"fast_pow({base}, {exp})"
        return f"std::pow({base}, {exp})"

    out = _POW_INT.sub(repl, expression)
    # Any remaining bare pow(...) calls that were not matched above get the namespace
    out = re.sub(r'(?<!std::)\bpow\s*\(', 'std::pow(', out)
    return out


def fix_exp(expr: str, funcs: Iterable[str] = FUNCS) -> str:
    """
    Add the std:: prefix to any math function calls that SymPy printed without it.
    Also collapses accidental double-namespace (std::std::) artifacts.
    """
    out = expr
    for fn in funcs:
        if fn in ("pow", "fast_pow"):
            continue  # handled separately by replace_pow
        out = make_pattern(fn).sub(f"std::{fn}(", out)

    out = re.sub(r'\bstd::\s*std::', 'std::', out)
    return out
