#pragma once
#include <array>
#include <cmath>
#include "../Functors/Functor.h"
#include "../Particle.h"

template <typename Particle_T>
class KryptonFunctorReference : public Functor<Particle_T> {
 public:
  explicit KryptonFunctorReference(
      double A, double a1, double a2, double a_m1, double b,
      double C6, double C8, double C10,
      bool newton3 = true,
      double cutoff = 0.0
  )
      : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b),
        _C6(C6), _C8(C8), _C10(C10),
        _newton3(newton3),
        _cutoff(cutoff)
  {
      // --- Dinamik recursion (makale / SymPy ile birebir) ---
      _C12 = _C6  * std::pow(_C10 / _C8, 3.0);
      _C14 = _C8  * std::pow(_C12 / _C10, 3.0);
      _C16 = _C10 * std::pow(_C14 / _C12, 3.0);
  }

  bool allowsNewton3() const  { return true; }
  bool usesNewton3()  const  { return _newton3; }

  void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
    const auto &r1 = p1.getR();
    const auto &r2 = p2.getR();

    const double dx = r1[0] - r2[0];
    const double dy = r1[1] - r2[1];
    const double dz = r1[2] - r2[2];

    double r2sq = dx*dx + dy*dy + dz*dz;
    constexpr double minR2 = 1e-24;
    if (r2sq < minR2) r2sq = minR2;

    
    if (_cutoff > 0.0) {
        const double cutoff2 = _cutoff * _cutoff;
        if (r2sq > cutoff2) return;
    }

    const double r     = std::sqrt(r2sq);
    const double invr  = 1.0 / r;
    const double invr2 = invr * invr;

    const double expAlphaTerm =
        std::exp(_a1 * r + _a2 * r*r + _a_m1 * invr);

    const double alphaTerm =
        _a1 + 2.0*_a2*r - _a_m1*invr2;

    const double firstTerm =
        (-_A * alphaTerm * expAlphaTerm) * invr;

    // ===== Tang–Toennies damping =====
    const double bdist = _b * r;
    const double expbr = std::exp(-bdist);

    // powers of b*r
    double bpow[17];
    bpow[0] = 1.0;
    for (int i = 1; i <= 16; ++i) bpow[i] = bpow[i-1] * bdist;

    constexpr double invF[17] = {
        1.,1.,1./2.,1./6.,1./24.,1./120.,1./720.,1./5040.,1./40320.,
        1./362880.,1./3628800.,1./39916800.,1./479001600.,
        1./6227020800.,1./87178291200.,1./1307674368000.,
        1./20922789888000.
    };

    double ksum[17];
    for (int i = 0; i <= 16; ++i)
        ksum[i] = bpow[i] * invF[i];

    const double ksumacc6  = ksum[0]+ksum[1]+ksum[2]+ksum[3]+ksum[4]+ksum[5]+ksum[6];
    const double ksumacc8  = ksumacc6  + ksum[7] + ksum[8];
    const double ksumacc10 = ksumacc8  + ksum[9] + ksum[10];
    const double ksumacc12 = ksumacc10 + ksum[11] + ksum[12];
    const double ksumacc14 = ksumacc12 + ksum[13] + ksum[14];
    const double ksumacc16 = ksumacc14 + ksum[15] + ksum[16];

    const double invr6  = invr2*invr2*invr2;
    const double invr8  = invr6 * invr2;
    const double invr10 = invr8 * invr2;
    const double invr12 = invr10 * invr2;
    const double invr14 = invr12 * invr2;
    const double invr16 = invr14 * invr2;

    const double term6  = _C6  * invr6  * (-6  + expbr*(6  *ksumacc6  + bdist*ksum[6]));
    const double term8  = _C8  * invr8  * (-8  + expbr*(8  *ksumacc8  + bdist*ksum[8]));
    const double term10 = _C10 * invr10 * (-10 + expbr*(10 *ksumacc10 + bdist*ksum[10]));
    const double term12 = _C12 * invr12 * (-12 + expbr*(12 *ksumacc12 + bdist*ksum[12]));
    const double term14 = _C14 * invr14 * (-14 + expbr*(14 *ksumacc14 + bdist*ksum[14]));
    const double term16 = _C16 * invr16 * (-16 + expbr*(16 *ksumacc16 + bdist*ksum[16]));

    const double secondTerm =
        (term6 + term8 + term10 + term12 + term14 + term16) * invr2;

    const double F_over_r = firstTerm + secondTerm;

    std::array<double,3> F{
        F_over_r * dx,
        F_over_r * dy,
        F_over_r * dz
    };

    p1.addF(F);
    if (_newton3) p2.subF(F);
  }

 private:
  double _A, _a1, _a2, _a_m1, _b;
  double _C6, _C8, _C10;
  double _C12, _C14, _C16;
  bool   _newton3;
  double _cutoff;
};
