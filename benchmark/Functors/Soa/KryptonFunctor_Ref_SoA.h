#pragma once
#include "../Particle.h"
#include <cstddef>
#include <cmath>
#include "FastPow.hpp"
template <class SoAView_T>
class KryptonFunctor_Ref_SoA {
public:
    using SoAView = SoAView_T;

    explicit KryptonFunctor_Ref_SoA(
        double A,
        double a1,
        double a2,
        double a_m1,
        double b,
        double C6,
        double C8,
        double C10,
        bool newton3 = true,
        double cutoff = 0.0
    )
        : _newton3(newton3)
        , _cutoff(cutoff)
        , _A(A)
        , _a1(a1)
        , _a2(a2)
        , _a_m1(a_m1)
        , _b(b)
        , _C6(C6)
        , _C8(C8)
        , _C10(C10)
    {}

    void SoAFunctor(SoAView &soa) {

        double *x  = soa.x.data();
        double *y  = soa.y.data();
        double *z  = soa.z.data();

        double *fx = soa.fx.data();
        double *fy = soa.fy.data();
        double *fz = soa.fz.data();

        const std::size_t N = soa.size();

        for (std::size_t i = 0; i < N; ++i) {
            double fix = 0.0, fiy = 0.0, fiz = 0.0;

            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];

            for (std::size_t j = i + 1; j < N; ++j) {

                const double dx = xi - x[j];
                const double dy = yi - y[j];
                const double dz = zi - z[j];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1e-24) continue;

                // -------- cutoff check --------
                if (_cutoff > 0.0) {
                    const double cutoff2 = _cutoff * _cutoff;
                    if (r2 > cutoff2) continue;
                }

                const double r     = std::sqrt(r2);
                const double inv_r = 1.0 / r;

                /* -----------------------------
                 * Potential terms
                 * ----------------------------- */

                // Short-range exponential part
                const double V_rep =
                    _A * std::exp(_a1*r + _a2*r2 + _a_m1*inv_r);

                // Tang–Toennies damping
                const double br = _b * r;
                const double exp_br = std::exp(-br);

                auto damp = [&](int n) {
                    double sum = 0.0;
                    double term = 1.0;
                    for (int k = 0; k <= n; ++k) {
                        if (k > 0) term *= br / k;
                        sum += term;
                    }
                    return 1.0 - exp_br * sum;
                };

                const double f6  = damp(6);
                const double f8  = damp(8);
                const double f10 = damp(10);

                const double r6  = std::pow(inv_r, 6);
                const double r8  = r6 * inv_r * inv_r;
                const double r10 = r8 * inv_r * inv_r;

                const double V_disp =
                    - f6  * _C6  * r6
                    - f8  * _C8  * r8
                    - f10 * _C10 * r10;

                /* -----------------------------
                 * Force = -dV/dr
                 * ----------------------------- */

                // dV_rep / dr
                const double dVrep =
                    V_rep * (_a1 + 2.0*_a2*r - _a_m1*inv_r*inv_r);

                // dV_disp / dr (numeric reference derivative)
                const double eps = 1e-6;
                const double rp  = r + eps;
                const double rm  = r - eps;

                auto V_disp_r = [&](double rr) {
                    const double ir = 1.0 / rr;
                    const double br_ = _b * rr;
                    const double exp_ = std::exp(-br_);

                    auto damp_ = [&](int n) {
                        double s = 0.0, t = 1.0;
                        for (int k = 0; k <= n; ++k) {
                            if (k > 0) t *= br_ / k;
                            s += t;
                        }
                        return 1.0 - exp_ * s;
                    };

                    return
                        - damp_(6)  * _C6  * std::pow(ir, 6)
                        - damp_(8)  * _C8  * std::pow(ir, 8)
                        - damp_(10) * _C10 * std::pow(ir,10);
                };

                const double dVdisp =
                    (V_disp_r(rp) - V_disp_r(rm)) / (2.0 * eps);

                const double Fmag = -(dVrep + dVdisp);

                const double Fx = Fmag * dx * inv_r;
                const double Fy = Fmag * dy * inv_r;
                const double Fz = Fmag * dz * inv_r;

                fix += Fx;
                fiy += Fy;
                fiz += Fz;

                if (_newton3) {
                    fx[j] -= Fx;
                    fy[j] -= Fy;
                    fz[j] -= Fz;
                }
            }

            fx[i] += fix;
            fy[i] += fiy;
            fz[i] += fiz;
        }
    }

private:
    bool _newton3;
    double _cutoff;

    double _A;
    double _a1;
    double _a2;
    double _a_m1;
    double _b;
    double _C6;
    double _C8;
    double _C10;
};
