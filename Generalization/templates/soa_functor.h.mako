#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max

template <class SoAView_T>
class ${classname}_SoA {
public:
    using SoAView = SoAView_T;

    explicit ${classname}_SoA(
% for name, value in parameters.items():
        double ${name},
% endfor
        bool newton3 = ${"true" if newton3_default else "false"}
    )
      : _newton3(newton3)
% for name, value in parameters.items():
      , _${name}(${name})
% endfor
% if classname.lower().startswith("mie"):
      , _C(0.0)
% endif
    {
% if classname.lower().startswith("mie"):
        // runtime C compute
        _C = (_n / (_n - _m)) * std::pow(_n / _m, _m / (_n - _m));
% endif
    }


    void SoAFunctor(SoAView &soa) {

        // --- RESTRICT pointers (SIMD-friendly) ---
        double * __restrict__ x  = soa.x.data();
        double * __restrict__ y  = soa.y.data();
        double * __restrict__ z  = soa.z.data();

        double * __restrict__ fx = soa.fx.data();
        double * __restrict__ fy = soa.fy.data();
        double * __restrict__ fz = soa.fz.data();

% if uses_mass:
        double * __restrict__ mass = soa.mass.data();
% endif

        const std::size_t N = soa.size();
        constexpr double EPS = ${eps_guard};

        // --- parameter aliases (loop-invariant) ---
% for name, value in parameters.items():
        const double ${name} = _${name};
% endfor
% if classname.lower().startswith("mie"):
        const double C = _C;
% endif

        for (std::size_t i = 0; i < N; ++i) {

            // --- local force accumulator for particle i ---
            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            // --- cache coordinates of i ---
            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];

% if uses_mass:
            const double p1m = mass[i];
% endif

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

% if uses_mass:
                const double p2m = mass[j];
% endif

% if temps_code:
                ${temps_code}
% endif

                const double Fmag = ${force_expr};

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
% for name, value in parameters.items():
    double _${name};
% endfor
% if classname.lower().startswith("mie"):
    double _C;
% endif
};
