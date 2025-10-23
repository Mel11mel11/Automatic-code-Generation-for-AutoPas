#pragma once
#ifndef GENERATED_MIE_HPP
#define GENERATED_MIE_HPP

#include <stdexcept>

namespace mie {

// declarations
double C(double n, double m);
double rmin(double sig, double n, double m);
double force_fast(double r, double eps, double sig, double n, double m);
double force_safe(double r, double eps, double sig, double n, double m);

// runtime switch — flag’siz kullanım 👇
inline double computeForce(double r, double eps, double sig, double n, double m, bool safeVersion = false) {
    if (safeVersion)
        return force_safe(r, eps, sig, n, m);
    else
        return force_fast(r, eps, sig, n, m);
}

} // namespace mie

#endif // GENERATED_MIE_HPP
