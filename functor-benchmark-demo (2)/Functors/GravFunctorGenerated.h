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
    explicit GravFunctorGenerated(double G, bool newton3=false) : _G(G), _newton3(newton3) {}

    void AoSFunctor(Particle_T& a, Particle_T& b) override {
        const auto& ra=a.getR(); const auto& rb=b.getR();
        const double dx=rb[0]-ra[0], dy=rb[1]-ra[1], dz=rb[2]-ra[2];

        double r2=dx*dx+dy*dy+dz*dz; if (r2 < 1e-24) r2 = 1e-24;
        const double r=std::sqrt(r2), invr=1.0/r;

        const double m1 = a.getMass(); // find mass of particle 1
        const double m2 = b.getMass(); // find mass of particle 2
        const double mag = grav::computeForce(r, _G, m1, m2); // |F|

        
        const std::array<double,3> F{ -mag*dx*invr, -mag*dy*invr, -mag*dz*invr };
        a.addF(F); if (_newton3) b.subF(F);
    }

private:
    double _G; bool _newton3;
};
