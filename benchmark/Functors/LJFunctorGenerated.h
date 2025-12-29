// LJFunctorGenerated.h
#pragma once

#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include "../Particle.h"

#include <array>
#include <cmath>

#include "../generated_files/generated_force.hpp"  // from python

template <typename Particle_T>
class LJFunctorGenerated : public Functor<Particle_T> {
public:
    explicit LJFunctorGenerated(
        double sigma,
        double epsilon,
        bool newton3 = true,
        double cutoff = 0.0
    )
        : _sigma(sigma),
          _epsilon(epsilon),
          _newton3(newton3),
          _cutoff(cutoff)
    {}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        const auto& r1 = p1.getR();
        const auto& r2 = p2.getR();

        const double dx = r1[0] - r2[0];
        const double dy = r1[1] - r2[1];
        const double dz = r1[2] - r2[2];

        double r2sq = dx*dx + dy*dy + dz*dz;

        constexpr double minR2 = 1e-8;
        if (r2sq < minR2) r2sq = minR2;

        // ---- cutoff (NO member cutoff2) ----
        if (_cutoff > 0.0) {
            const double cutoff2 = _cutoff * _cutoff;
            if (r2sq > cutoff2) return;
        }

        const double r = std::sqrt(r2sq);
        const double inv_r = 1.0 / r;

        const double mag = lj::computeForce(r, _epsilon, _sigma);

        const std::array<double,3> F{
            mag * dx * inv_r,
            mag * dy * inv_r,
            mag * dz * inv_r
        };

        p1.addF(F);
        if (_newton3) p2.subF(F);
    }

    bool usesNewton3() const { return _newton3; }
    bool allowsNewton3() const  { return true; }

private:
    double _sigma;
    double _epsilon;
    bool   _newton3;
    double _cutoff;
};
