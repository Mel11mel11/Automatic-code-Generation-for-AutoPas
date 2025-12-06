#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>

template <class SoAView_T>
class GravityFunctor_Gen_SoA {
public:
    using SoAView = SoAView_T;

    explicit GravityFunctor_Gen_SoA(
        double G,
        bool newton3 = true
    )
        : _newton3(newton3)
        , _G(G)
    {
    }

   
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
            for (std::size_t j = _newton3 ? i + 1 : 0; j < N; ++j) {

                if (i == j) continue;

                double dx = x[j] - x[i];
                double dy = y[j] - y[i];
                double dz = z[j] - z[i];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1e-24) r2 = 1e-24;

                double r = std::sqrt(r2);
                double inv_r = 1.0 / r;

                double p1m = mass[i];
                double p2m = mass[j];

                const double G = _G;


                double Fmag = -G*fast_pow(inv_r, 2)*p1m*p2m;

                double Fx = Fmag * dx * inv_r;
                double Fy = Fmag * dy * inv_r;
                double Fz = Fmag * dz * inv_r;

                fx[i] += Fx;
                fy[i] += Fy;
                fz[i] += Fz;

                if (_newton3) {
                    fx[j] -= Fx;
                    fy[j] -= Fy;
                    fz[j] -= Fz;
                }
            }
        }

    } // SoAFunctor

private:
    bool _newton3;
    double _G;
};
