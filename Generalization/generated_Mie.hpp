#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>

template <class Particle_T>
class MieFunctor_Gen : public Functor<Particle_T> {
public:
  MieFunctor_Gen(double sigma, double epsilon, int n, int m, double C, bool newton3=true) : _newton3(newton3), _sigma2(sigma*sigma), _n(n), _m(m), _C(C), _epsilon(epsilon) {}

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
    double sig_pow_n = 1.0, sig_pow_m = 1.0;
    { // σ^2 tabanı
      double s2 = _sigma2;
      int hn = _n/2, hm = _m/2;
      for(int i=0;i<hn;i++) sig_pow_n *= s2;
      for(int i=0;i<hm;i++) sig_pow_m *= s2;
    }
    double inv_r_pow_n2 = 1.0, inv_r_pow_m2 = 1.0;
    {{
      int pn = (_n+2)/2, pm = (_m+2)/2;
      for(int i=0;i<pn;i++) inv_r_pow_n2 *= inv_r2;
      for(int i=0;i<pm;i++) inv_r_pow_m2 *= inv_r2;
    }}
    const double term_n = _n * sig_pow_n * inv_r_pow_n2;
    const double term_m = _m * sig_pow_m * inv_r_pow_m2;
    const double fac = _C * _epsilon * (term_n - term_m) * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;
    a.addF(Fx, Fy, Fz);
    if (_newton3) b.subF(Fx, Fy, Fz);
  }

private:
  bool _newton3;
double _sigma2 = 0.0;
  int _n = 0;
  int _m = 0;
  double _C = 0.0;
  double _epsilon = 0.0;
};
