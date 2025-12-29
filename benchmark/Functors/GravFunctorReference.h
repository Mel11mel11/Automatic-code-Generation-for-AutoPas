#pragma once
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include <cmath>

// Reference implementation of the Newtonian gravitational force
template <typename Particle_T>
class GravFunctorReference : public Functor<Particle_T> {
public:
    explicit GravFunctorReference(
        double gravConst,
        bool newton3 = false,
        double cutoff = 0.0
    )
        : _gravConst(gravConst),
          _newton3(newton3),
          _cutoff(cutoff) {}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        using namespace arrayMath::literals;

        const auto& ra = p1.getR();
        const auto& rb = p2.getR();
        auto dr = ra - rb;  // displacement vector

        double r2 = arrayMath::dot(dr, dr);
        constexpr double minR2 = 1e-24;
        if (r2 < minR2) r2 = minR2;

        
        if (_cutoff > 0.0) {
            const double cutoff2 = _cutoff * _cutoff;
            if (r2 > cutoff2) return;
        }

        const double r = std::sqrt(r2);
        const double invr3 = 1.0 / (r * r2);

        const double m1 = p1.getMass();
        const double m2 = p2.getMass();

        const double mag = -_gravConst * m1 * m2 * invr3;
        auto F = dr * mag;

        p1.addF(F);
        if (_newton3) p2.subF(F);
    }

    bool usesNewton3() const  { return _newton3; }
    bool allowsNewton3() const  { return true; }

private:
    double _gravConst;
    bool   _newton3;
    double _cutoff;
};
