#include <cmath>
#include "generated_force.hpp"

namespace lj {
inline double computeForceMag_noPow_impl(double r, double eps24, double sigma6) {
    // Guard against r=0 from caller side if needed; this routine assumes r>0.
    double inv_r  = 1.0 / r;
    double inv_r2 = inv_r * inv_r;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;

    double s6_over_r6 = sigma6 * inv_r6;                // (sigma/r)^6
    return eps24 * (2.0 * s6_over_r6 * s6_over_r6 - s6_over_r6) * inv_r;
}

double computeForceMag_noPow(double r, double eps24, double sigma6) {
    // Optional tiny clamp for stability (aynı minR2 mantığını r tabanında uygulamak istersen):
    // constexpr double minR = 1e-4;  // örnek (r^2=1e-8 ile uyumlu)
    // if (r < minR) r = minR;

    return computeForceMag_noPow_impl(r, eps24, sigma6);
}
} // namespace lj
