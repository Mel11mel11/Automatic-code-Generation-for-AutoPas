#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max

template <class SoAView_T>
class GravityFunctor_Gen_O100_SoA {
public:
    using SoAView = SoAView_T;

    explicit GravityFunctor_Gen_O100_SoA(
        double G,
        bool newton3 = true,
        double cutoff = 0.0
    )
      : _newton3(newton3)
      , _cutoff(cutoff)
      , _G(G)
    {
    }

    void SoAFunctor(SoAView &soa) {

        // --- SoA raw pointers ---
        double * __restrict__ x  = soa.x.data();
        double * __restrict__ y  = soa.y.data();
        double * __restrict__ z  = soa.z.data();

        double * __restrict__ fx = soa.fx.data();
        double * __restrict__ fy = soa.fy.data();
        double * __restrict__ fz = soa.fz.data();

        double * __restrict__ mass = soa.mass.data();

        const std::size_t N = soa.size();
        constexpr double EPS = 1e-24;

        const double cutoff  = _cutoff;
        const double cutoff2 = cutoff * cutoff;

        const double G = _G;

        for (std::size_t i = 0; i < N; ++i) {

            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];

            const double p1m = mass[i];

            if (_newton3) {
                // =====================================================
                // Newton3 = ON  → half-shell (i < j), symmetric update
                // =====================================================
                for (std::size_t j = i + 1; j < N; ++j) {

                    const double dx = xi - x[j];
                    const double dy = yi - y[j];
                    const double dz = zi - z[j];

                    double r2 = dx*dx + dy*dy + dz*dz;
                    r2 = std::max(r2, EPS);

                    if (cutoff > 0.0 && r2 > cutoff2) continue;

                    const double r = std::sqrt(r2);
                    const double inv_r = 1.0 / r;

                    const double p2m = mass[j];


                    const double Fmag = -G*fast_pow(inv_r, 2)*p1m*p2m;

                    const double Fx = Fmag * dx * inv_r;
                    const double Fy = Fmag * dy * inv_r;
                    const double Fz = Fmag * dz * inv_r;

                    fix += Fx;
                    fiy += Fy;
                    fiz += Fz;

                    fx[j] -= Fx;
                    fy[j] -= Fy;
                    fz[j] -= Fz;
                }

            } else {
                // =====================================================
                // Newton3 = OFF → full-shell (i != j), single-sided
                // =====================================================
                for (std::size_t j = 0; j < N; ++j) {
                    if (j == i) continue;

                    const double dx = xi - x[j];
                    const double dy = yi - y[j];
                    const double dz = zi - z[j];

                    double r2 = dx*dx + dy*dy + dz*dz;
                    r2 = std::max(r2, EPS);

                    if (cutoff > 0.0 && r2 > cutoff2) continue;

                    const double r = std::sqrt(r2);
                    const double inv_r = 1.0 / r;

                    const double p2m = mass[j];


                    const double Fmag = -G*fast_pow(inv_r, 2)*p1m*p2m;

                    fix += Fmag * dx * inv_r;
                    fiy += Fmag * dy * inv_r;
                    fiz += Fmag * dz * inv_r;
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

    double _G;
};
