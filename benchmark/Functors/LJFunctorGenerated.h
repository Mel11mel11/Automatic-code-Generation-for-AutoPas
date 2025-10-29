// LJFunctorGenerated.h
#pragma once
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include "Functor.h"
#include "../Particle.h"
#include <array>
#include <cmath>
#include "../generated_files/generated_force.hpp"   // header from my python script

template <typename Particle_T>
class LJFunctorGenerated : public Functor<Particle_T> {
public:
    LJFunctorGenerated(double sigma, double epsilon, bool newton3=false)
      : _sigma(sigma), _epsilon(epsilon), _newton3(newton3) {}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        const auto& r1 = p1.getR();
        const auto& r2 = p2.getR();
        const double dx = r1[0]-r2[0]; // abs function?
        const double dy = r1[1]-r2[1];
        const double dz = r1[2]-r2[2];

        double r2sq = dx*dx + dy*dy + dz*dz;
        if (r2sq < 1e-24) r2sq = 1e-24;      // softening (stability)
        const double r = std::sqrt(r2sq);
        //const double mag = lj::computeForce(r, _epsilon, _sigma);
        //const std::array<double,3> F{ -mag*dx, -mag*dy, -mag*dz };  // <-- eksi eklendi
        const double mag = lj::computeForce(r, _epsilon, _sigma);
        const double inv_r = 1.0 / r;
        const std::array<double,3> F{ mag * dx * inv_r,
                              mag * dy * inv_r,
                              mag * dz * inv_r };
        //const double inv_r = 1.0 / r;
        //const std::array<double,3> F{ mag*dx*inv_r, mag*dy*inv_r, mag*dz*inv_r }; // F*r_12/r

        p1.addF(F);
        if (_newton3) p2.subF(F);
    }
    
bool usesNewton3() const { return _newton3; }
void setNewton3(bool state) { _newton3 = state; }

private:
    double _sigma, _epsilon;
    bool _newton3;
};
