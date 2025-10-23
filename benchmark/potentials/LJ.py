import sympy as sp
from sympy.printing.c import ccode

# Step 1- Define variables 
r, epsilon, sigma = sp.symbols('r epsilon sigma')
# Step 2- Define Lennard-Jones potential 
V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)
# Step 3- Compute the force 
F = -sp.diff(V, r)
F_simplified = sp.simplify(F) # Simplify the expression
print(f"\033[1;94mSimplified Force F(r) = {F_simplified}\033[0m\n")
# Step 4- Generate C++ code 
cpp_expr = ccode(F_simplified) # only knows <cmath> so need regulations for  std::pow or smilar etc.
# math functions for -c++ version
cpp_expr = (cpp_expr
            .replace("pow(",  "std::pow(")
            .replace("sqrt(", "std::sqrt(")
            .replace("fabs(", "std::fabs(")
            .replace("cbrt(", "std::cbrt(")
            .replace("exp(", "std::exp("))
         
# Step 5- Compose C++ function 
cpp_function = f"""
#include <cmath> // for math functions
#include "generated_force.hpp" // include the header file 

// original auto-generated Lennard-Jones force from sympy python script
namespace lj {{
double computeForce(double r, double epsilon, double sigma) {{
    return {cpp_expr};
}}
}}
"""
#  Step 6- Header file 
header_code = """#pragma once
namespace lj {
    double computeForce(double r, double epsilon, double sigma);
}
"""

# Step 7- Write files
with open("../generated_files/generated_force.cpp", "w") as f:
    f.write(cpp_function)

with open("../generated_files/generated_force.hpp", "w") as f:
    f.write(header_code)
print("\033[1;33mC++ code generated successfully 🟢 "
      "Look into generated_force.cpp and .hpp files!\033[0m\n")
