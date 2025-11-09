#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>

template <class Particle_T>
class GravityFunctor_Gen : public Functor<Particle_T> {
public:
  GravityFunctor_Gen(double G, double m1, double m2, bool newton3=true) : _newton3(newton3), _G(G), _m1(m1), _m2(m2) {}

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

    const double inv_r2 = 1.0 / r2;
    const double inv_r  = 1.0 / std::sqrt(r2);
    const double inv_r3 = inv_r2 * inv_r;
    const double C = _G * _m1 * _m2;
    const double fac = -C * inv_r3 * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;
    a.addF(Fx, Fy, Fz);
    if (_newton3) b.subF(Fx, Fy, Fz);
  }

private:
  bool _newton3;
double _G = 0.0;
  double _m1 = 0.0;
  double _m2 = 0.0;
};
