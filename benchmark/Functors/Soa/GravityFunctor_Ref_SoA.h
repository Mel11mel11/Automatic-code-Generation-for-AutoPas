#pragma once
#include "../Particle.h"
#include <cstddef>
#include <cmath>

template <class SoAView_T>
class GravityFunctor_Ref_SoA {
public:
    using SoAView = SoAView_T;

    explicit GravityFunctor_Ref_SoA(
        double G,
        bool newton3 = true
    )
        : _newton3(newton3)
        , _G(G)
    {}

    void SoAFunctor(SoAView &soa) {
        double *x  = soa.x.data();
        double *y  = soa.y.data();
        double *z  = soa.z.data();

        double *fx = soa.fx.data();
        double *fy = soa.fy.data();
        double *fz = soa.fz.data();

        double *mass = soa.mass.data();

        const std::size_t N = soa.size();

        for (std::size_t i = 0; i < N; ++i) {
            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];
            const double mi = mass[i];

            for (std::size_t j = i + 1; j < N; ++j) {
                const double dx = xi - x[j];
                const double dy = yi - y[j];
                const double dz = zi - z[j];

                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1e-24) continue;

                const double inv_r = 1.0 / std::sqrt(r2);
                const double inv_r3 = inv_r * inv_r * inv_r;

                // Attractive gravity:
                // F_i = -G * mi * mj * (r_i - r_j) / |r|^3
                const double fac = -_G * mi * mass[j] * inv_r3;

                const double fxij = dx * fac;
                const double fyij = dy * fac;
                const double fzij = dz * fac;

                fix += fxij;
                fiy += fyij;
                fiz += fzij;

                if (_newton3) {
                    fx[j] -= fxij;
                    fy[j] -= fyij;
                    fz[j] -= fzij;
                }
            }

            fx[i] += fix;
            fy[i] += fiy;
            fz[i] += fiz;
        }
    }

private:
    bool _newton3;
    double _G;
};
