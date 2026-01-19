#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max

template <class SoAView_T>
class TestFunctor_Gen_Opt00000_SoA {
public:
    using SoAView = SoAView_T;

    explicit TestFunctor_Gen_Opt00000_SoA(
        bool newton3 = true,
        double cutoff = 0.0
    )
      : _newton3(newton3)
      , _cutoff(cutoff)
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



                // Contract: Fmag == fr == -dU/dr (scalar)
                const double Fmag = -2/inv_r;

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

};
