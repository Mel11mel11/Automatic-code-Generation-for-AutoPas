#pragma once
#include <cmath>

namespace lj {
inline double computeForce(double r, double epsilon, double sigma) {
    if (r < 1e-12) r = 1e-12;            // for security
    const double sr = sigma / r;         // (σ/r)
    const double sr2 = sr*sr;
    const double sr6 = sr2*sr2*sr2;
    const double sr12 = sr6*sr6;
    return 24.0*epsilon * (2.0*sr12 - sr6) / r;
}
}