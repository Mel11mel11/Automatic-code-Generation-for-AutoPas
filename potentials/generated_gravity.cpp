#include <cmath>
#include "generated_gravity.hpp"

// Auto-generated Gravitational Force
namespace grav {
double computeForce(double r, double G, double m1, double m2) {
    return -G*m1*m2/std::pow(r, 2);
}
}
