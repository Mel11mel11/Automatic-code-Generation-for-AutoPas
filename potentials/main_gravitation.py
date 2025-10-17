import sympy as sp
from sympy.printing.c import ccode

# === Step 1: Define variables ===
r, G, m1, m2 = sp.symbols('r G m1 m2')

# === Step 2: Define gravitational potential ===
V = -G * m1 * m2 / r
print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)

# === Step 3: Compute force (negative derivative) ===
F = -sp.diff(V, r)
F_simplified = sp.simplify(F)
print(f"\033[1;94mSimplified Force F(r) = {F_simplified}\033[0m\n")
# === Step 4: Generate C++ code ===
cpp_expr = ccode(F_simplified)

# C++ tarafında ad-qualify: pow → std::pow, sqrt → std::sqrt, fabs → std::fabs
cpp_expr = (cpp_expr
            .replace("pow(",  "std::pow(")
            .replace("sqrt(", "std::sqrt(")
            .replace("fabs(", "std::fabs("))

# === Step 5: Compose C++ source ===
cpp_function = f"""#include <cmath>
#include "generated_gravity.hpp"

// Auto-generated Gravitational Force
namespace grav {{
double computeForce(double r, double G, double m1, double m2) {{
    return {cpp_expr};
}}
}}
"""

# === Step 6: Compose header ===
cpp_header = """#pragma once
namespace grav {
    double computeForce(double r, double G, double m1, double m2);
}
"""

# === Step 7: Write files ===
with open("generated_gravity.cpp", "w") as f:
    f.write(cpp_function)

with open("generated_gravity.hpp", "w") as f:
    f.write(cpp_header)

print("\033[1;33mC++ code generated successfully 🟢 "
      "Look into generated_gravity.cpp and .hpp files!\033[0m\n")

