import sympy as sp
from  sympy import symbols, diff
from sympy.printing.c import ccode
from sympy.printing import print_ccode
import os

# define the variables one by one
r, epsilon, sigma = symbols('r epsilon sigma')

# define Lennard-Jones potential 
V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) # U(r) = 4\epsilon[(\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^6]

print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)

# compute the force with negative differentiation 
F = - diff(V, r)
# print_ccode(F)  # print C code version of the expression
F_simplified = sp.simplify(F) # for time complexity reduction simplifying the expression
# print_ccode(F_simplified)  # print C code version of the expression
print(f"\033[1;94mSimplified Force F(r) = {F_simplified}\033[0m\n") # see the output 

# generate C++ code with ccode() function
cpp_expr = ccode(F_simplified) # only knows <cmath> so need regulations for  std::pow or smilar etc.

# math functions for c++ version
cpp_expr = (cpp_expr
            .replace("pow(",  "std::pow(")
            .replace("sqrt(", "std::sqrt(")
            .replace("fabs(", "std::fabs(")
            .replace("cbrt(", "std::cbrt(")
            .replace("exp(", "std::exp("))
         
# compose C++ function for cpp document

cpp_function = f"""
#include <cmath> // for math functions
#include "generated_force.hpp" // include the header file 

// auto-generated Lennard-Jones force from sympy python script
namespace lj {{
double computeForce(double r, double epsilon, double sigma) {{
    return {cpp_expr};
}}
}}
"""
# header file .hpp content
header_code = f"""#pragma once // for include guard
namespace lj {{
    double computeForce(double r, double epsilon, double sigma);
}}
"""

os.makedirs("../generated_files", exist_ok=True)
# write files
with open("../generated_files/generated_force.cpp", "w") as d:
    d.write(cpp_function)

with open("../generated_files/generated_force.hpp", "w") as d:
    d.write(header_code)

print("\033[1;33mC++ code generated successfully 🟢 "
      "Look into generated_force.cpp and .hpp files!\033[0m\n")
