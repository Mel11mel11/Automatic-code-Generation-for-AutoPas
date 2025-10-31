
import sympy as sp
from sympy.printing.c import ccode

# Step 1- Define variables 
A,a_1,R,a_2,a_11,b,k,C , sigma = sp.symbols('r epsilon sigma')
# Step 2- Define Lennard-Jones potential 
# V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
print("\033[1;36mPotential V(r) = %s\033[0m\n" % V)
# Step 3- Compute the force 
F = -sp.diff(V, r)
F_simplified = sp.simplify(F)
