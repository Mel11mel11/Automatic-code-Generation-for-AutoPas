#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>

template <class SoAView_T>
class LJFunctor_Gen_SoA {
public:
    using SoAView = SoAView_T;

    explicit LJFunctor_Gen_SoA(
        double sigma,
        double epsilon,
        bool newton3 = true
    )
      : _newton3(newton3)
      , _sigma(sigma)
      , _epsilon(epsilon)
    {
    }

    // =============================
    //        SOA FUNCTOR
    // =============================
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
            // ALWAYS j = i+1 (AoS ile birebir)
            for (std::size_t j = i + 1; j < N; ++j) {

                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];

                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 < 1e-24) r2 = 1e-24;

                double r = std::sqrt(r2);
                double inv_r = 1.0 / r;

                double p1m = mass[i];
                double p2m = mass[j];

                const double sigma = _sigma;
                const double epsilon = _epsilon;

const double x0 = fast_pow(sigma, 6);

                double Fmag = -24*epsilon*fast_pow(inv_r, 13)*x0*(fast_pow(r, 6) - 2*x0);

                double Fx = Fmag * dx * inv_r;
                double Fy = Fmag * dy * inv_r;
                double Fz = Fmag * dz * inv_r;

                fx[i] += Fx; fy[i] += Fy; fz[i] += Fz;

                if (_newton3) {
                    fx[j] -= Fx;
                    fy[j] -= Fy;
                    fz[j] -= Fz;
                }
            }
        }
    }

private:
    bool _newton3;
    double _sigma;
    double _epsilon;
};
