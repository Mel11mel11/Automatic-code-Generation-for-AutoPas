#pragma once
#include <array>
#include <cmath>
#include "../Functors/Functor.h"
#include "../Particle.h"

template <typename Particle_T>
class KryptonFunctorReference : public Functor<Particle_T> {
 public:
  explicit KryptonFunctorReference(double A, double a1, double a2, double a_m1, double b,
                                   double C6, double C8, double C10,
                                   bool newton3 = true)
      : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b),
        _C6(C6), _C8(C8), _C10(C10),
        _newton3(newton3)
  {
      // --- Dinamik recursion (SymPy / makale ile birebir) ---
      _C12 = _C6  * std::pow(_C10 / _C8, 3.0);
      _C14 = _C8  * std::pow(_C12 / _C10, 3.0);
      _C16 = _C10 * std::pow(_C14 / _C12, 3.0);
  }

  bool allowsNewton3() const { return true; }
  bool usesNewton3()  const { return _newton3; }

  void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
    const auto &r1 = p1.getR();
    const auto &r_2 = p2.getR();

    const double dx = r1[0] - r_2[0];
    const double dy = r1[1] - r_2[1];
    const double dz = r1[2] - r_2[2];

    double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 < 1e-24) r2 = 1e-24;

    const double r = std::sqrt(r2);
    const double invr = 1.0 / r;
    const double invr2 = invr * invr;

   
    const double expAlphaTerm = std::exp(_a1 * r + _a2 * r*r + _a_m1 * invr);
    const double alphaTerm = _a1 + 2.0*_a2*r - _a_m1*invr2;

    const double firstTerm = (-_A * alphaTerm * expAlphaTerm) * invr;

    const double bdist  = _b * r;
    const double expbr  = std::exp(-bdist); // buraya bak krypton için mutlaka belki alggemeine vereifacxhung bulursun


    const double b2 = bdist * bdist;
    const double b3 = b2 * bdist;
    const double b4 = b3 * bdist;
    const double b5 = b4 * bdist;
    const double b6 = b5 * bdist;
    const double b7 = b6 * bdist;
    const double b8 = b7 * bdist;
    const double b9 = b8 * bdist;
    const double b10 = b9 * bdist;
    const double b11 = b10 * bdist;
    const double b12 = b11 * bdist;
    const double b13 = b12 * bdist;
    const double b14 = b13 * bdist;
    const double b15 = b14 * bdist;
    const double b16 = b15 * bdist;


    constexpr double invF[17] = {
        1.,1.,1./2.,1./6.,1./24.,1./120.,1./720.,1./5040.,1./40320.,
        1./362880.,1./3628800.,1./39916800.,1./479001600.,1./6227020800.,
        1./87178291200.,1./1307674368000.,1./20922789888000.
    };


    const double ksum0  = 1.0;
    const double ksum1  = bdist;
    const double ksum2  = b2  * invF[2];
    const double ksum3  = b3  * invF[3];
    const double ksum4  = b4  * invF[4];
    const double ksum5  = b5  * invF[5];
    const double ksum6  = b6  * invF[6];
    const double ksum7  = b7  * invF[7];
    const double ksum8  = b8  * invF[8];
    const double ksum9  = b9  * invF[9];
    const double ksum10 = b10 * invF[10];
    const double ksum11 = b11 * invF[11];
    const double ksum12 = b12 * invF[12];
    const double ksum13 = b13 * invF[13];
    const double ksum14 = b14 * invF[14];
    const double ksum15 = b15 * invF[15];
    const double ksum16 = b16 * invF[16];

    const double ksumacc6  = ksum0+ksum1+ksum2+ksum3+ksum4+ksum5+ksum6;
    const double ksumacc8  = ksumacc6 + ksum7 + ksum8;
    const double ksumacc10 = ksumacc8 + ksum9 + ksum10;
    const double ksumacc12 = ksumacc10 + ksum11 + ksum12;
    const double ksumacc14 = ksumacc12 + ksum13 + ksum14;
    const double ksumacc16 = ksumacc14 + ksum15 + ksum16;


    const double invr6  = invr2*invr2*invr2;
    const double invr8  = invr6 * invr2;
    const double invr10 = invr8 * invr2;
    const double invr12 = invr10 * invr2;
    const double invr14 = invr12 * invr2;
    const double invr16 = invr14 * invr2;


    const double term6  = _C6  * invr6  * (-6  + expbr*(6  *ksumacc6  + bdist*ksum6));
    const double term8  = _C8  * invr8  * (-8  + expbr*(8  *ksumacc8  + bdist*ksum8));
    const double term10 = _C10 * invr10 * (-10 + expbr*(10 *ksumacc10 + bdist*ksum10));
    const double term12 = _C12 * invr12 * (-12 + expbr*(12 *ksumacc12 + bdist*ksum12));
    const double term14 = _C14 * invr14 * (-14 + expbr*(14 *ksumacc14 + bdist*ksum14));
    const double term16 = _C16 * invr16 * (-16 + expbr*(16 *ksumacc16 + bdist*ksum16));

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16) * invr2;

    const double F_over_r = firstTerm + secondTerm;

    const double Fx = F_over_r * dx;
    const double Fy = F_over_r * dy;
    const double Fz = F_over_r * dz;

    std::array<double,3> F{Fx, Fy, Fz};

    p1.addF(F);
    if (_newton3)
      p2.subF(F);
  }

 private:
  double _A, _a1, _a2, _a_m1, _b;
  double _C6, _C8, _C10;
  double _C12, _C14, _C16;
  bool _newton3;
};
