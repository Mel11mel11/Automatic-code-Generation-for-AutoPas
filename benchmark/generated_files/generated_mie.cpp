#include <cmath>
#include "generated_mie.hpp"

namespace mie {

double computeForce(double r, double eps, double sig, double n, double m) {
    if (r < 1e-12) r = 1e-12;
    return eps*std::pow(m, m/(m - n))*std::pow(n, -n/(m - n))*std::pow(r, -m - n - 2)*(m*std::pow(r, n + 1)*std::pow(sig, m) - n*std::pow(r, m + 1)*std::pow(sig, n))/(m - n);
}

} // namespace mie
