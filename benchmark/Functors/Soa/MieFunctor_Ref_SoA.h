#pragma once
#include "../Particle.h"
#include <cstddef>
#include <cmath>
#include <cassert>

template <class SoAView_T>
class MieFunctor_Ref_SoA {
public:
    using SoAView = SoAView_T;

    explicit MieFunctor_Ref_SoA(
        double sigma,
        double epsilon,
        int n,
        int m,
        bool newton3 = true
    )
        : _newton3(newton3)
        , _sigma(sigma)
        , _sigma2(sigma * sigma)
        , _epsilon(epsilon)
        , _n(n)
        , _m(m)
    {
        // Basic sanity checks (optional but helpful)
        assert(_n > 0 && _m > 0);
        assert(_n > _m);

        // Normalization constant C
        _C = (_n / double(_n - _m)) *
             std::pow(double(_n) / double(_m),
                      double(_m) / double(_n - _m));
    }

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

                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1e-24) continue;

                const double inv_r2 = 1.0 / r2;
                const double s2_over_r2 = _sigma2 * inv_r2;

                // We need (sigma/r)^n and (sigma/r)^m for integer n,m (odd allowed).
                // Compute baseEvenPow = (sigma^2/r^2)^(floor(exp/2)) via multiplications.
                // If exp is odd, multiply once more by (sigma/r) = sigma * (1/sqrt(r2)).
                const double inv_r = 1.0 / std::sqrt(r2);
                const double sr = _sigma * inv_r;  // (sigma/r)

                const double sn = pow_sigma_over_r_int(sr, s2_over_r2, _n);
                const double sm = pow_sigma_over_r_int(sr, s2_over_r2, _m);

                // Force factor for vector form:
                // F_vec = - C*eps * ( n*(sigma/r)^n - m*(sigma/r)^m ) * r_vec / r^2
                const double fac = -_C * _epsilon * (_n * sn - _m * sm) * inv_r2;

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
    // Computes (sigma/r)^exp for integer exp using:
    // (sigma/r)^exp = (sigma^2/r^2)^(exp/2) * (sigma/r)^(exp%2)
    static inline double pow_sigma_over_r_int(double sr, double s2_over_r2, int exp) {
        const int half = exp / 2;

        double evenPart = 1.0;
        for (int k = 0; k < half; ++k) evenPart *= s2_over_r2;

        if (exp & 1) {
            return evenPart * sr; // multiply once more by (sigma/r)
        }
        return evenPart;
    }

    bool _newton3;
    double _sigma;
    double _sigma2;
    double _epsilon;
    double _C;
    int _n;
    int _m;
};
