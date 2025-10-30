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
   explicit GravFunctorGenerated(double gravConst, bool newton3=false)
 : _gravConst(gravConst), _newton3(newton3) {}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
    const auto& r1 = p1.getR();
    const auto& r2 = p2.getR();

    const double dx = r1[0] - r2[0];
    const double dy = r1[1] - r2[1];
    const double dz = r1[2] - r2[2];

    double r2sq = dx*dx + dy*dy + dz*dz;
    if (r2sq < 1e-24) r2sq = 1e-24;
    const double r = std::sqrt(r2sq);

    // computeForce = -G*m1*m2 / r^2  → vektöre çevirmek için /r
    const double coeff = grav::computeForce(r, _gravConst, p1.getMass(), p2.getMass()) / r;

    std::array<double,3> F = { coeff * dx, coeff * dy, coeff * dz };

    p1.addF(F);
    if (_newton3) p2.subF(F);  // 
}

bool usesNewton3() const  { return _newton3; }

private:
    double _gravConst ; bool _newton3;
};
