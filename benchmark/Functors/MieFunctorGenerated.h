#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <array>
#include <cmath>
#include "../generated_files/generated_mie.hpp"

template <typename Particle_T>
class MieFunctorGenerated : public Functor<Particle_T> {
public:
    MieFunctorGenerated(
        double sigma,
        double epsilon,
        int n,
        int m,
        bool newton3 = true,
        double cutoff = 0.0
    )
        : _sigma(sigma),
          _epsilon(epsilon),
          _n(n),
          _m(m),
          _newton3(newton3),
          _cutoff(cutoff)
    {}

    void AoSFunctor(Particle_T& a, Particle_T& b) override {
        const auto& ra = a.getR();
        const auto& rb = b.getR();

        const double dx = ra[0] - rb[0];
        const double dy = ra[1] - rb[1];
        const double dz = ra[2] - rb[2];

        double r2 = dx*dx + dy*dy + dz*dz;
        constexpr double minR2 = 1e-24;
        if (r2 < minR2) r2 = minR2;

        // ---- cutoff (AutoPas-style, local r² check) ----
        if (_cutoff > 0.0) {
            const double cutoff2 = _cutoff * _cutoff;
            if (r2 > cutoff2) return;
        }

        const double r    = std::sqrt(r2);
        const double invr = 1.0 / r;

        const double mag =
            mie::computeForce(r, _epsilon, _sigma, _n, _m);

        const std::array<double,3> F{
            mag * dx * invr,
            mag * dy * invr,
            mag * dz * invr
        };

        a.addF(F);
        if (_newton3) b.subF(F);
    }

    bool usesNewton3() const { return _newton3; }
    bool allowsNewton3() const { return true; }

private:
    double _sigma;
    double _epsilon;
    int    _n;
    int    _m;
    bool   _newton3;
    double _cutoff;
};
