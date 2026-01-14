#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max
#include "FastPow.hpp"
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


                const double tt_x = b*r;
        const double exp_m_tt_x = std::exp(-tt_x);
        double tt_arr[17];
        tt_arr[0] = 1.0;
        double tt_pow = 1.0;
        double tt_inv_fact = 1.0;
        double tt_sum = 1.0;
        for (int k = 1; k <= 16; ++k) {
          tt_pow *= tt_x;
          tt_inv_fact /= double(k);
          tt_sum += tt_pow * tt_inv_fact;
          tt_arr[k] = tt_sum;
        }
        const double tt5 = tt_arr[5];
        const double tt6 = tt_arr[6];
        const double tt7 = tt_arr[7];
        const double tt8 = tt_arr[8];
        const double tt9 = tt_arr[9];
        const double tt10 = tt_arr[10];
        const double tt11 = tt_arr[11];
        const double tt12 = tt_arr[12];
        const double tt13 = tt_arr[13];
        const double tt14 = tt_arr[14];
        const double tt15 = tt_arr[15];
        const double tt16 = tt_arr[16];
        const double x0 = std::exp(-b*r);
        const double x1 = b*x0;
        const double x2 = C6*x1;
        const double x3 = fast_pow(C10, 3)/fast_pow(C8, 3);
        const double x4 = fast_pow(C10, 6)*fast_pow(C6, 3)/fast_pow(C8, 8);
        const double x5 = fast_pow(C10, 10)*fast_pow(C6, 6)/fast_pow(C8, 15);
        const double x6 = fast_pow(r, 2);

                // Contract: Fmag == fr == -dU/dr (scalar)
                const double Fmag = -A*(a1 + 2*a2*r - a_m1/x6)*std::exp(a1*r + a2*x6 + a_m1/r) + C10*x1*(tt10 - tt9)/fast_pow(r, 10) + 10*C10*(tt10*x0 - 1)/fast_pow(r, 11) + 6*C6*(tt6*x0 - 1)/fast_pow(r, 7) + 12*C6*x3*(tt12*x0 - 1)/fast_pow(r, 13) - C8*x1*(tt7 - tt8)/fast_pow(r, 8) + 8*C8*(tt8*x0 - 1)/fast_pow(r, 9) - x2*(tt5 - tt6)/fast_pow(r, 6) - x2*x3*(tt11 - tt12)/fast_pow(r, 12) - x1*x4*(tt13 - tt14)/fast_pow(r, 14) + 14*x4*(tt14*x0 - 1)/fast_pow(r, 15) - x1*x5*(tt15 - tt16)/fast_pow(r, 16) + 16*x5*(tt16*x0 - 1)/fast_pow(r, 17);

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
