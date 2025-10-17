// generated_mie.cpp (auto-produced by mie_codegen.py)
#include <cmath>
#include <stdexcept>
#include "generated_mie.hpp"

namespace mie {

static inline void validate(double r, double eps, double sig, double n, double m) {
    if (!(r > 0.0))  throw std::invalid_argument("Mie: r must be > 0");
    if (!(eps > 0.0)) throw std::invalid_argument("Mie: epsilon must be > 0");
    if (!(sig > 0.0)) throw std::invalid_argument("Mie: sigma must be > 0");
    if (!(n > 0.0 && m > 0.0)) throw std::invalid_argument("Mie: n,m must be > 0");
    if (!(n > m)) throw std::invalid_argument("Mie: require n > m for a proper well");
}

double C(double n, double m) {
    // C(n,m) = n*pow(n/m, m/(-m + n))/(-m + n)
    return n*pow(n/m, m/(-m + n))/(-m + n);
}

double rmin(double sig, double n, double m) {
    // r_min = sig*pow(n/m, 1.0/(-m + n))
    return sig*pow(n/m, 1.0/(-m + n));
}

double force_fast(double r, double eps, double sig, double n, double m) {
    // Common subexpressions:
    const double x0 = 1.0/(m - n);
    // Force:
    return eps*pow(m, m*x0)*pow(n, -n*x0)*pow(r, -m - n - 2)*x0*(m*pow(r, n + 1)*pow(sig, m) - n*pow(r, m + 1)*pow(sig, n));
}

double force_safe(double r, double eps, double sig, double n, double m) {
    validate(r, eps, sig, n, m);
    const double x0 = 1.0/(m - n);
    return eps*pow(m, m*x0)*pow(n, -n*x0)*pow(r, -m - n - 2)*x0*(m*pow(r, n + 1)*pow(sig, m) - n*pow(r, m + 1)*pow(sig, n));
}

} // namespace mie
