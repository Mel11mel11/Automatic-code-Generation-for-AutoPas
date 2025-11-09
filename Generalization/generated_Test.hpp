#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>

template <class Particle_T>
class TestFunctor_Gen : public Functor<Particle_T> {
public:
  TestFunctor_Gen(bool newton3=true) : _newton3(newton3) {}

  bool allowsNewton3() const { return true; }
  bool usesNewton3()  const { return _newton3; }

  inline void AoSFunctor(Particle_T& a, Particle_T& b) override {
    const auto& ra = a.getR(); const auto& rb = b.getR();
    const double dx = ra[0]-rb[0];
    const double dy = ra[1]-rb[1];
    const double dz = ra[2]-rb[2];
    const double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 == 0.0) return;

    const double mask = 1.0;

    const double fac = -2.0 * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;
    a.addF(Fx, Fy, Fz);
    if (_newton3) b.subF(Fx, Fy, Fz);
  }

private:
  bool _newton3;

};
