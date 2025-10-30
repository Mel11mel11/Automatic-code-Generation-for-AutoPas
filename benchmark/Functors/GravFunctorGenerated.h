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

    void AoSFunctor(Particle_T& a, Particle_T& b) override {
        
        const auto& ra=a.getR(); const auto& rb=b.getR();
        const double dx=ra[0]-rb[0], dy=ra[1]-rb[1], dz=ra[2]-rb[2];

        double r2=dx*dx+dy*dy+dz*dz; if (r2 < 1e-24) r2 = 1e-24;
        const double r=std::sqrt(r2), invr=1.0/r;

        const double m1 = a.getMass();
        const double m2 = b.getMass(); // 
        const double mag = grav::computeForce(r, _gravConst, m1, m2);
        // ---- DEBUG print (scientific) ----
{
    //auto oldflags = std::cout.flags();
  //  auto oldprec  = std::cout.precision();

  /*  std::cout.setf(std::ios::scientific);
    std::cout.precision(12);
    std::cout << "[TEST] r=" << r
              << "  G="   << _gravConst
              << "  m1="  << m1
              << "  m2="  << m2
              << "  mag=" << mag
              << std::endl;

    std::cout.flags(oldflags);
    std::cout.precision(oldprec);*/
}
// ---- end DEBUG ----
 // |F|
        const std::array<double,3> F{ mag*dx*invr, mag*dy*invr, mag*dz*invr };

        a.addF(F); 
        if (_newton3) b.subF(F);
    }

private:
    double _gravConst ; bool _newton3;
};
