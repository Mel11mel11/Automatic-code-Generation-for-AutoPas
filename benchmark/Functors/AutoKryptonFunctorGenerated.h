#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <typename Particle_T>
class AutoKryptonFunctorGenerated : public Functor<Particle_T> {
public:
    explicit AutoKryptonFunctorGenerated(double A, double a1, double a2, double a_m1, double bo,
                                         double C6, double C8, double C10,
                                         bool newton3)
        : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _bo(bo),
          _C6(C6), _C8(C8), _C10(C10),
          _newton3(newton3) {}

    void AoSFunctor(Particle_T& a, Particle_T& b) override {

        using namespace arrayMath::literals;
        const auto& ra = a.getR();
        const auto& rb = b.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        // r^2 and r; guard against r -> 0
        constexpr double EPS = 1e-24;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;
        const double r     = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        // Parameter aliases
        const double A     = _A;
        const double a1    = _a1;
        const double a2    = _a2;
        const double a_m1  = _a_m1;
        const double bo    = _bo;
        const double C6    = _C6;
        const double C8    = _C8;
        const double C10   = _C10;

        // --- helper sums for exp series ---
        auto truncated_exp = [](double z, int N) {
            double term = 1.0;
            double sum  = 1.0; // k=0 term
            for (int k = 1; k <= N; ++k) {
                term *= z / static_cast<double>(k);
                sum  += term;
            }
            return sum;
        };

        auto truncated_exp_weighted = [](double z, int N) {
            double term = 1.0;
            double sum  = 0.0; // k=0 term is 0 anyway
            for (int k = 1; k <= N; ++k) {
                term *= z / static_cast<double>(k);
                sum  += static_cast<double>(k) * term;
            }
            return sum;
        };

        const double z   = bo * r;

        const double S6  = truncated_exp(z, 6);
        const double S6w = truncated_exp_weighted(z, 6);
        const double S8  = truncated_exp(z, 8);
        const double S8w = truncated_exp_weighted(z, 8);
        const double S10  = truncated_exp(z, 10);
        const double S10w = truncated_exp_weighted(z, 10);
        const double S12  = truncated_exp(z, 12);
        const double S12w = truncated_exp_weighted(z, 12);
        const double S14  = truncated_exp(z, 14);
        const double S14w = truncated_exp_weighted(z, 14);
        const double S16  = truncated_exp(z, 16);
        const double S16w = truncated_exp_weighted(z, 16);

        // --- codegen: force magnitude F(r, params) ---
        const double Fmag =
            -A * (a1 + 2.0 * a2 * r - a_m1 / fast_pow(r, 2))
              * std::exp(a1 * r + a2 * fast_pow(r, 2) + a_m1 / r)

            + fast_pow(C10, 10) * fast_pow(C6, 6)
              * ( bo * std::exp(-bo * r) * S16
                  - std::exp(-bo * r) * (S16w / r) )
              / ( std::pow(C8, 15) * std::pow(r, 16) )

            - 16.0 * fast_pow(C10, 10) * fast_pow(C6, 6)
              * ( 1.0 - std::exp(-bo * r) * S16 )
              / ( std::pow(C8, 15) * std::pow(r, 17) )

            + fast_pow(C10, 6) * fast_pow(C6, 3)
              * ( bo * std::exp(-bo * r) * S14
                  - std::exp(-bo * r) * (S14w / r) )
              / ( fast_pow(C8, 8) * std::pow(r, 14) )

            - 14.0 * fast_pow(C10, 6) * fast_pow(C6, 3)
              * ( 1.0 - std::exp(-bo * r) * S14 )
              / ( fast_pow(C8, 8) * std::pow(r, 15) )

            + fast_pow(C10, 3) * C6
              * ( bo * std::exp(-bo * r) * S12
                  - std::exp(-bo * r) * (S12w / r) )
              / ( fast_pow(C8, 3) * fast_pow(r, 12) )

            - 12.0 * fast_pow(C10, 3) * C6
              * ( 1.0 - std::exp(-bo * r) * S12 )
              / ( fast_pow(C8, 3) * fast_pow(r, 13) )

            + C10
              * ( bo * std::exp(-bo * r) * S10
                  - std::exp(-bo * r) * (S10w / r) )
              / fast_pow(r, 10)

            - 10.0 * C10
              * ( 1.0 - std::exp(-bo * r) * S10 )
              / fast_pow(r, 11)

            + C6
              * ( bo * std::exp(-bo * r) * S6
                  - std::exp(-bo * r) * (S6w / r) )
              / fast_pow(r, 6)

            - 6.0 * C6
              * ( 1.0 - std::exp(-bo * r) * S6 )
              / fast_pow(r, 7)

            + C8
              * ( bo * std::exp(-bo * r) * S8
                  - std::exp(-bo * r) * (S8w / r) )
              / fast_pow(r, 8)

            - 8.0 * C8
              * ( 1.0 - std::exp(-bo * r) * S8 )
              / fast_pow(r, 9);

        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double, 3> F{fx, fy, fz};

        a.addF(F);
        if (_newton3) {
            b.subF(F);
        }
    }

    bool allowsNewton3() const  { return true; }
    bool usesNewton3()   const  { return _newton3; }

private:
    double _A;
    double _a1;
    double _a2;
    double _a_m1;
    double _bo;
    double _C6;
    double _C8;
    double _C10;
    bool   _newton3;
};
