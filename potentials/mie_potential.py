# mie_codegen.py
import sympy as sp
from sympy.printing.c import ccode

n, m, r, eps, sig = sp.symbols('n m r eps sig', positive=True)

# Normalization: C(n,m) = n/(n-m) * (n/m)^(m/(n-m))
C = (n/(n - m)) * (n/m)**(m/(n - m))

sr = sig / r
V = C * eps * (sr**n - sr**m)          # Mie potential
F = -sp.diff(V, r) 
print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)                    
F = sp.simplify(F)
print(f"\033[1;94mSimplified Force F(r) = {F}\033[0m\n")

# Also emit helpers
rmin_expr = sig * (n/m)**(1/(n - m))

#  Step 2: CSE to reduce pow calls 
repl, reduced = sp.cse([F], optimizations='basic')
F_reduced = reduced[0]

# Prepare code pieces
repl_code = [f"const double {ccode(sym)} = {ccode(expr)};"
             for (sym, expr) in repl]
F_code    = ccode(F_reduced)
C_code    = ccode(C)
rmin_code = ccode(rmin_expr)

# Step 3: Write header 
hpp = f"""#pragma once
#ifndef GENERATED_MIE_HPP
#define GENERATED_MIE_HPP

#include <stdexcept>

namespace mie {{
double C(double n, double m);
double rmin(double sig, double n, double m);
double force_fast(double r, double eps, double sig, double n, double m);
double force_safe(double r, double eps, double sig, double n, double m);
}} // namespace mie

#endif // GENERATED_MIE_HPP
"""

with open("generated_mie.hpp", "w") as f:
    f.write(hpp)

# Step 4: Write implementation 
cpp = f"""// generated_mie.cpp (auto-produced by mie_codegen.py)
#include <cmath>
#include <stdexcept>
#include "generated_mie.hpp"

namespace mie {{

static inline void validate(double r, double eps, double sig, double n, double m) {{
    if (!(r > 0.0))  throw std::invalid_argument("Mie: r must be > 0");
    if (!(eps > 0.0)) throw std::invalid_argument("Mie: epsilon must be > 0");
    if (!(sig > 0.0)) throw std::invalid_argument("Mie: sigma must be > 0");
    if (!(n > 0.0 && m > 0.0)) throw std::invalid_argument("Mie: n,m must be > 0");
    if (!(n > m)) throw std::invalid_argument("Mie: require n > m for a proper well");
}}

double C(double n, double m) {{
    // C(n,m) = {C_code}
    return {C_code};
}}

double rmin(double sig, double n, double m) {{
    // r_min = {rmin_code}
    return {rmin_code};
}}

double force_fast(double r, double eps, double sig, double n, double m) {{
    // Common subexpressions:
{"\n".join("    " + line for line in repl_code)}
    // Force:
    return {F_code};
}}

double force_safe(double r, double eps, double sig, double n, double m) {{
    validate(r, eps, sig, n, m);
{"\n".join("    " + line for line in repl_code)}
    return {F_code};
}}

}} // namespace mie
"""

with open("generated_mie.cpp", "w") as f:
    f.write(cpp)

print("\033[1;33mC++ code generated successfully 🟢 "
      "Look into generated_mie.cpp and .hpp files!\033[0m\n")
