
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <array>
#include <cmath>
#include "../generated_files/generated_mie.hpp"

template <typename Particle_T>
class MieFunctorGenerated : public Functor<Particle_T> {
public:
    MieFunctorGenerated(double sigma, double epsilon, int n, int m, bool newton3)
    : _sigma(sigma), _epsilon(epsilon), _n(n), _m(m), _newton3(newton3) {}

    void AoSFunctor(Particle_T& a, Particle_T& b) override {
        const auto& ra=a.getR();
        const auto& rb=b.getR();
        const double dx=ra[0]-rb[0], dy=ra[1]-rb[1], dz=ra[2]-rb[2];
        double r2=dx*dx+dy*dy+dz*dz; 

        if(r2<1e-24) r2=1e-24;
        
        const double r=std::sqrt(r2), invr=1.0/r;
        const double mag = mie::computeForce(r, _epsilon, _sigma, _n, _m);
        const std::array<double,3> F_safe{mag*dx*invr, mag*dy*invr, mag*dz*invr};

       
        a.addF(F_safe);
        if (_newton3) b.subF(F_safe);
    }
    bool usesNewton3() const  { return _newton3; }
    bool allowsNewton3() const { return true; }

private: double _sigma,_epsilon; int _n,_m; bool _newton3;
};
