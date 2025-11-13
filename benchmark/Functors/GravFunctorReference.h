#pragma once
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include <cmath>
#include <iostream>

// Reference implementation of the Newtonian gravitational force
template <typename Particle_T>
class GravFunctorReference : public Functor<Particle_T> {
public:
explicit GravFunctorReference(double gravConst, bool newton3)
    : _gravConst(gravConst), _newton3(newton3) {}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        using namespace arrayMath::literals;  // enables array math ops (+, -, *)

        const auto& ra = p1.getR();
        const auto& rb = p2.getR();
        auto dr = ra - rb;  // displacement vector

        double r2 = arrayMath::dot(dr, dr);
        if (r2 < 1e-24) r2 = 1e-24;  // softening for numerical stability
        double r = std::sqrt(r2);
        double invr3 = 1.0 / (r2 * r);

        double m1 = p1.getMass();
        double m2 = p2.getMass();
        double mag =- _gravConst * m1 * m2 * invr3;
        auto F = dr * mag;  // scale vector
        p1.addF(F);
        if (_newton3) p2.subF(F);
    }
    bool usesNewton3() const { return _newton3; }
    bool allowsNewton3() const { return true; }

private:
    double  _gravConst;
    bool _newton3;
};
