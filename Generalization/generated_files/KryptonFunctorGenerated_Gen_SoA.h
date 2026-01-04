#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max

template <class SoAView_T>
class KryptonFunctorGenerated_Gen_SoA {
public:
    using SoAView = SoAView_T;

    explicit KryptonFunctorGenerated_Gen_SoA(
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
    {
    }

    void SoAFunctor(SoAView &soa) {

        // --- RESTRICT pointers (SIMD-friendly) ---
        double * __restrict__ x  = soa.x.data();
        double * __restrict__ y  = soa.y.data();
        double * __restrict__ z  = soa.z.data();

        double * __restrict__ fx = soa.fx.data();
        double * __restrict__ fy = soa.fy.data();
        double * __restrict__ fz = soa.fz.data();


        const std::size_t N = soa.size();
        constexpr double EPS = 1e-24;

        // --- loop-invariant config ---
        const double cutoff = _cutoff;
        const double cutoff2 = cutoff * cutoff;

        // --- parameter aliases (loop-invariant) ---
        const double A = _A;
        const double a1 = _a1;
        const double a2 = _a2;
        const double a_m1 = _a_m1;
        const double b = _b;
        const double C6 = _C6;
        const double C8 = _C8;
        const double C10 = _C10;

        for (std::size_t i = 0; i < N; ++i) {

            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];


            for (std::size_t j = i + 1; j < N; ++j) {

                const double dx = xi - x[j];
                const double dy = yi - y[j];
                const double dz = zi - z[j];

                double r2 = dx * dx + dy * dy + dz * dz;
                r2 = std::max(r2, EPS);

                if (cutoff > 0.0 && r2 > cutoff2) continue;

                // --- define r and inv_r consistently ---
                const double r = std::sqrt(r2);
                const double inv_r = 1.0 / r;


                const double x0 = fast_pow(b, 11);
        const double x1 = b*r;
        const double x2 = std::exp(-x1);
        const double x3 = fast_pow(b, 7);
        const double x4 = C6*x2;
        const double x5 = fast_pow(b, 9);
        const double x6 = fast_pow(b, 13);
        const double x7 = fast_pow(C10, 3)/fast_pow(C8, 3);
        const double x8 = fast_pow(b, 15);
        const double x9 = fast_pow(C10, 6)*fast_pow(C6, 3)/fast_pow(C8, 8);
        const double x10 = fast_pow(C10, 10)*fast_pow(C6, 6)/fast_pow(C8, 15);
        const double x11 = fast_pow(r, 2);
        const double x12 = fast_pow(r, 7);
        const double x13 = fast_pow(b, 6)*fast_pow(r, 6);
        const double x14 = fast_pow(b, 2)*x11;
        const double x15 = fast_pow(b, 3)*fast_pow(r, 3);
        const double x16 = fast_pow(b, 4)*fast_pow(r, 4);
        const double x17 = fast_pow(b, 5)*fast_pow(r, 5);
        const double x18 = fast_pow(r, 9);
        const double x19 = fast_pow(b, 8)*fast_pow(r, 8);
        const double x20 = x12*x3;
        const double x21 = fast_pow(r, 11);
        const double x22 = fast_pow(b, 10)*fast_pow(r, 10);
        const double x23 = x18*x5;
        const double x24 = fast_pow(r, 13);
        const double x25 = fast_pow(b, 12)*fast_pow(r, 12);
        const double x26 = x0*x21;
        const double x27 = fast_pow(r, 15);
        const double x28 = fast_pow(b, 14)*fast_pow(r, 14);
        const double x29 = x24*x6;

                // Contract: Fmag == fr == -dU/dr (scalar)
                const double Fmag = -A*(a1 + 2*a2*r - a_m1/x11)*std::exp(a1*r + a2*x11 + a_m1/r) + (1.0/3628800.0)*C10*x0*x2 + (1.0/362880.0)*C10*(x2*(3628800*x1 + 5040*x13 + 1814400*x14 + 604800*x15 + 151200*x16 + 30240*x17 + 90*x19 + 720*x20 + x22 + 10*x23 + 3628800) - 3628800)/x21 + (1.0/39916800.0)*C6*x7*(x2*(479001600*x1 + 665280*x13 + 239500800*x14 + 79833600*x15 + 19958400*x16 + 3991680*x17 + 11880*x19 + 95040*x20 + 132*x22 + 1320*x23 + x25 + 12*x26 + 479001600) - 479001600)/x24 + (1.0/120.0)*C6*(x2*(720*x1 + x13 + 360*x14 + 120*x15 + 30*x16 + 6*x17 + 720) - 720)/x12 + (1.0/40320.0)*C8*x2*x5 + (1.0/5040.0)*C8*(x2*(40320*x1 + 56*x13 + 20160*x14 + 6720*x15 + 1680*x16 + 336*x17 + x19 + 8*x20 + 40320) - 40320)/x18 + (1.0/20922789888000.0)*fast_pow(b, 17)*x10*x2 + (1.0/87178291200.0)*x2*x8*x9 + (1.0/720.0)*x3*x4 + (1.0/479001600.0)*x4*x6*x7 + (1.0/6227020800.0)*x9*(x2*(87178291200*x1 + 121080960*x13 + 43589145600*x14 + 14529715200*x15 + 3632428800*x16 + 726485760*x17 + 2162160*x19 + 17297280*x20 + 24024*x22 + 240240*x23 + 182*x25 + 2184*x26 + x28 + 14*x29 + 87178291200) - 87178291200)/x27 + (1.0/1307674368000.0)*x10*(x2*(fast_pow(b, 16)*fast_pow(r, 16) + 20922789888000*x1 + 29059430400*x13 + 10461394944000*x14 + 3487131648000*x15 + 871782912000*x16 + 174356582400*x17 + 518918400*x19 + 4151347200*x20 + 5765760*x22 + 57657600*x23 + 43680*x25 + 524160*x26 + 16*x27*x8 + 240*x28 + 3360*x29 + 20922789888000) - 20922789888000)/fast_pow(r, 17);

                // Convert to vector: F_vec = fr * r_vec / r
                const double Fx = Fmag * dx * inv_r;
                const double Fy = Fmag * dy * inv_r;
                const double Fz = Fmag * dz * inv_r;

                fix += Fx;
                fiy += Fy;
                fiz += Fz;

                // Newton3: only apply to j if enabled
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
