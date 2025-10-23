#include <cmath>
#include <stdexcept>
#include "generated_mie.hpp"

namespace mie {

static void validate(double r, double eps, double sig, double n, double m) {
    if (!(r > 0.0)) throw std::invalid_argument("Mie: r must be > 0");
    if (!(eps > 0.0)) throw std::invalid_argument("Mie: epsilon must be > 0");
    if (!(sig > 0.0)) throw std::invalid_argument("Mie: sigma must be > 0");
    if (!(n > 0.0 && m > 0.0)) throw std::invalid_argument("Mie: n,m must be > 0");
    if (!(n > m)) throw std::invalid_argument("Mie: require n > m for a proper well");
}

double C(double n, double m) {
    return n*std::pow(n/m, m/(-m + n))/(-m + n);
}

double rmin(double sig, double n, double m) {
    return sig*std::pow(n/m, 1.0/(-m + n));
}

double force_fast(double r, double eps, double sig, double n, double m) {
    if (r < 1e-12) r = 1e-12;
    const double x0 = 1.0/(m - n);
    return eps*std::pow(m, m*x0)*std::pow(n, -n*x0)*std::pow(r, -m - n - 2)*x0*
           (m*std::pow(r, n + 1)*std::pow(sig, m) - n*std::pow(r, m + 1)*std::pow(sig, n));
}

double force_safe(double r, double eps, double sig, double n, double m) {
    validate(r, eps, sig, n, m);
    if (r < 1e-12) r = 1e-12;
    const double x0 = 1.0/(m - n);
    return eps*std::pow(m, m*x0)*std::pow(n, -n*x0)*std::pow(r, -m - n - 2)*x0*
           (m*std::pow(r, n + 1)*std::pow(sig, m) - n*std::pow(r, m + 1)*std::pow(sig, n));
}

} // namespace mie
