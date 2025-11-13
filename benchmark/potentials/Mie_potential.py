import sympy as sp
from sympy.printing.c import ccode
import os
# define symbolic variables and assume each of them is positive (greater than 0)
n, m, r, eps, sig = sp.symbols('n m r eps sig', positive=True)

# Normalization constant C(n,m)
C = (n/(n - m)) * (n/m)**(m/(n - m))

# potential V(r)
sr = sig / r
V  = C * eps * (sr**n - sr**m)
print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)

# force F(r) = -dV/dr
F = -sp.diff(V, r)
F = sp.simplify(F)
print(f"\033[1;94mSimplified Force F(r) = {F}\033[0m\n")

F_code = ccode(F)

F_code = (F_code
            .replace("pow(",  "std::pow(")
            .replace("sqrt(", "std::sqrt(")
            .replace("fabs(", "std::fabs(")
            .replace("cbrt(", "std::cbrt(")
            .replace("exp(", "std::exp("))
         
hpp = f"""#pragma once
#include <cmath>

namespace mie {{

// auto-generated Mie potential force from sympy python script 
double computeForce(double r, double eps, double sig, double n, double m);

}} // namespace mie
"""
os.makedirs("../generated_files", exist_ok=True)

with open("../generated_files/generated_mie.hpp", "w") as d:
    d.write(hpp)

cpp = f"""#include <cmath>
#include "generated_mie.hpp"

namespace mie {{

double computeForce(double r, double eps, double sig, double n, double m) {{
    if (r < 1e-12) r = 1e-12;
    return {F_code};
}}

}}
"""

with open("../generated_files/generated_mie.cpp", "w") as d:
    d.write(cpp)

print("\033[1;33mC++ code generated successfully 🟢 "
      "Look into generated_mie.cpp and .hpp files!\033[0m\n")
