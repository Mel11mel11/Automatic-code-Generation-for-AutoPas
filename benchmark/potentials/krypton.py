
import sympy as sp
from sympy.printing.c import ccode
import os
# This formulas is based on the 8.Eq. from https://www.hsu-hh.de/thermodynamik/wp-content/uploads/sites/741/2019/06/State-of-the-art-ab-initio-potential-energy-curve-for-the-krypton-atom-pair-and-thermophysical-properties-of-dilute-krypton-gas.pdf
# Page:10

R = sp.Symbol('R', positive=True)

A, a1, a2, a_m1, b = sp.symbols('A a1 a2 a_m1 b', real=True)
C6, C8, C10, C12, C14, C16 = sp.symbols('C6 C8 C10 C12 C14 C16', real=True)

V_exp = A * sp.exp(a1*R + a2*R**2 + a_m1/R)

C_list = [C6, C8, C10, C12, C14, C16]
n_values = list(range(3, 9))  # 3 to 8

sum_disp = 0
for Cn, n in zip(C_list, n_values):
    inner = sum((b*R)**k / sp.factorial(k) for k in range(0, 2*n + 1))
    expr  = (Cn / R**(2*n)) * (1 - sp.exp(-b*R) * inner)
    sum_disp += expr

V = V_exp - sum_disp

print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)
#sp.pretty_print(V)

F = -sp.diff(V, R)
F_simplified = sp.simplify(F)

print(f"\033[1;94mSimplified Force F(r) = {F_simplified}\033[0m\n") # see the output
#sp.pretty_print(F_simplified)

F_code = ccode(F_simplified)  # uses exp, pow, etc.

# map C math to <cmath> with std::
F_code = (F_code
            .replace("pow(",  "std::pow(")
            .replace("sqrt(", "std::sqrt(")
            .replace("fabs(", "std::fabs(")
            .replace("cbrt(", "std::cbrt(")
            .replace("exp(",  "std::exp("))


cpp_function = f"""
#include <cmath>   
#include "generated_krypton_force.hpp"

// Auto-generated Krypton force (Eq.8) from SymPy
namespace krypton {{
double computeForce(double R,
                          double A, double a1, double a2, double a_m1, double b,
                          double C6, double C8, double C10, double C12, double C14, double C16) {{
    return {F_code};
}}
}} 
"""

header_code = """#pragma once
namespace krypton {
    double computeForce(double R,
                          double A, double a1, double a2, double a_m1, double b,
                          double C6, double C8, double C10, double C12, double C14, double C16);
}
"""


os.makedirs("../generated_files", exist_ok=True)

with open("../generated_files/generated_krypton_force.cpp", "w") as f:
    f.write(cpp_function)

with open("../generated_files/generated_krypton_force.hpp", "w") as f:
    f.write(header_code)

print("\n\033[1;33mC++ code generated successfully 🟢 "
      "See generated_krypton_force.cpp / .hpp\033[0m\n")
