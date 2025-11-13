import sympy as sp
from sympy.printing.c import ccode
import os
# define variables ===
r, G, m1, m2 = sp.symbols('r G m1 m2')

# the formula for gravitational potential V(r) = -G * m1 * m2 / r
V = -G * m1 * m2 / r
print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)

# compute the force
F = -sp.diff(V, r)
F_simplified = sp.simplify(F)
print(f"\033[1;94mSimplified Force F(r) = {F_simplified}\033[0m\n")
# generate C++ code from the simplified force expression
F_code = ccode(F_simplified)

# replace the math functions to use std:: from <cmath>
F_code = (F_code
            .replace("pow(",  "std::pow(")
            .replace("sqrt(", "std::sqrt(")
            .replace("fabs(", "std::fabs("))

# compose C++ source
cpp_function = f"""#include <cmath>
#include "generated_gravity.hpp"

// Auto-generated Gravitational Force
namespace grav {{
double computeForce(double r, double G, double m1, double m2) {{
    return {F_code};
}}
}}
"""

# Compose header file
cpp_header = """#pragma once
namespace grav {
    double computeForce(double r, double G, double m1, double m2);
}
"""

# writing process for .hpp and .cpp
with open("../generated_files/generated_gravity.cpp", "w") as f:
    f.write(cpp_function)

with open("../generated_files/generated_gravity.hpp", "w") as f:
    f.write(cpp_header)

print("\033[1;33mC++ code generated successfully 🟢 "
      "Look into generated_gravity.cpp and .hpp files!\033[0m\n")

