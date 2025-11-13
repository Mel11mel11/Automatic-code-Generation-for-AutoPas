
#include <cmath> // for math functions
#include "generated_force.hpp" // include the header file 

// auto-generated Lennard-Jones force from sympy python script
namespace lj {
double computeForce(double r, double epsilon, double sigma) {
    return 24*epsilon*std::pow(sigma, 6)*(-std::pow(r, 6) + 2*std::pow(sigma, 6))/std::pow(r, 13);
}
}
