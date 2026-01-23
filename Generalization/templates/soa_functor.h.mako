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
        bool newton3 = ${"true" if newton3_default else "false"},
        double cutoff = 0.0
    )
      : _newton3(newton3)
      , _cutoff(cutoff)
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

% if use_mass:
        double * __restrict__ mass = soa.mass.data();
% endif

        const std::size_t N = soa.size();
        constexpr double EPS = ${eps_guard};

        // --- loop-invariant config ---
        const double cutoff = _cutoff;
        const double cutoff2 = cutoff * cutoff;

        // --- parameter aliases (loop-invariant) ---
% for name, value in parameters.items():
        const double ${name} = _${name};
% endfor
% if classname.lower().startswith("mie"):
        const double C = _C;
% endif

        for (std::size_t i = 0; i < N; ++i) {

            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];

% if use_mass:
            const double p1m = mass[i];
% endif

            for (std::size_t j = i + 1; j < N; ++j) {

                const double dx = xi - x[j];
                const double dy = yi - y[j];
                const double dz = zi - z[j];

                double r2 = dx * dx + dy * dy + dz * dz;
                r2 = std::max(r2, EPS);

% if cutoff_enabled:
                if (cutoff > 0.0 && r2 > cutoff2) continue;
% endif

                // --- define r and inv_r consistently ---
                const double r = std::sqrt(r2);
                const double inv_r = 1.0 / r;

% if use_mass:
                const double p2m = mass[j];
% endif

% if temps_code:
                ${temps_code}
% endif

                // Contract: Fmag == fr == -dU/dr (scalar)
                const double Fmag = ${force_expr};

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

% for name, value in parameters.items():
    double _${name};
% endfor
% if classname.lower().startswith("mie"):
    double _C;
% endif
};
