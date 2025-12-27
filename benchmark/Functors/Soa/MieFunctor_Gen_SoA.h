#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max

template <class SoAView_T>
class MieFunctor_Gen_SoA {
public:
    using SoAView = SoAView_T;

    explicit MieFunctor_Gen_SoA(
        double sigma,
        double epsilon,
        double n,
        double m,
        bool newton3 = true
    )
      : _newton3(newton3)
      , _sigma(sigma)
      , _epsilon(epsilon)
      , _n(n)
      , _m(m)
      , _C(0.0)
    {
        // runtime C compute
        _C = (_n / (_n - _m)) * std::pow(_n / _m, _m / (_n - _m));
    }

    // =============================
    //        SOA FUNCTOR
    // =============================
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

        // --- parameter aliases (loop-invariant) ---
        const double sigma = _sigma;
        const double epsilon = _epsilon;
        const double n = _n;
        const double m = _m;
        const double C = _C;

        for (std::size_t i = 0; i < N; ++i) {

            // --- local force accumulator for particle i ---
            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            // --- cache coordinates of i ---
            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];


            // Newton3 true -> j starts at i+1
            // Newton3 false -> we still start at i+1 but we must update both i and j
            for (std::size_t j = i + 1; j < N; ++j) {

                // cache coordinates of j (helps memory subsystem a bit)
                const double xj = x[j];
                const double yj = y[j];
                const double zj = z[j];

                const double dx = xi - xj;
                const double dy = yi - yj;
                const double dz = zi - zj;

                double r2 = dx * dx + dy * dy + dz * dz;
                // branchless guard (SIMD-friendly)
                r2 = std::max(r2, EPS);

                // cheaper: compute inv_r directly; avoid separate r variable
                const double inv_r = 1.0 / std::sqrt(r2);


                const double x0 = inv_r*sigma;

                const double Fmag = -C*epsilon*inv_r*(m*std::pow(x0, m) - n*std::pow(x0, n));

                // Keep this form to match your symbolic expression usage (dx * inv_r)
                const double Fx = Fmag * dx * inv_r;
                const double Fy = Fmag * dy * inv_r;
                const double Fz = Fmag * dz * inv_r;

                // accumulate force on i locally
                fix += Fx;
                fiy += Fy;
                fiz += Fz;

                // update force on j:
                // - Newton3: subtract (pair computed once)
                // - no Newton3: we still need to add the opposite contribution to j
                fx[j] -= Fx;
                fy[j] -= Fy;
                fz[j] -= Fz;
            }

            // write back force for i ONCE
            fx[i] += fix;
            fy[i] += fiy;
            fz[i] += fiz;
        }
    }

private:
    bool _newton3;
    double _sigma;
    double _epsilon;
    double _n;
    double _m;
    double _C;
};
