// Functors/GravFunctorGenerated.h
#pragma once
#include "Functor.h"
#include "../Particle.h"
#include <array>
#include <cmath>
#include "../generated_files/generated_gravity.hpp"

template <typename Particle_T>
class GravFunctorGenerated : public Functor<Particle_T> {
public:
    explicit GravFunctorGenerated(
        double gravConst,
        bool newton3 = false,
        double cutoff = 0.0
    )
        : _gravConst(gravConst),
          _newton3(newton3),
          _cutoff(cutoff) {}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        const auto& r1 = p1.getR();
        const auto& r2 = p2.getR();

        const double dx = r1[0] - r2[0];
        const double dy = r1[1] - r2[1];
        const double dz = r1[2] - r2[2];

        double r2sq = dx*dx + dy*dy + dz*dz;
        constexpr double minR2 = 1e-24;
        if (r2sq < minR2) r2sq = minR2;

        // ---- cutoff (AutoPas-style, local) ----
        if (_cutoff > 0.0) {
            const double cutoff2 = _cutoff * _cutoff;
            if (r2sq > cutoff2) return;
        }

        const double r = std::sqrt(r2sq);

        // computeForce = -G*m1*m2 / r^2
        const double m1 = p1.getMass();
        const double m2 = p2.getMass();

        // grav::computeForce(r, G, m1, m2) = -G*m1*m2 / r^2
        // we multiply by (dx/r) later → divide by r once more
        const double coeff =
            grav::computeForce(r, _gravConst, m1, m2) / r;

        std::array<double,3> F{
            coeff * dx,
            coeff * dy,
            coeff * dz
        };

        p1.addF(F);
        if (_newton3) p2.subF(F);
    }

    bool usesNewton3() const { return _newton3; }
    bool allowsNewton3() const  { return true; }

private:
    double _gravConst;
    bool   _newton3;
    double _cutoff;
};
