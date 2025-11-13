#pragma once
#include <cstddef>
#include <utility>
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include "../Functors/Functor.h"
#include <array>
#include <cmath>
#include "../Particle.h"

// Credits: https://github.com/AutoPas/AutoPas/blob/feat/3xa/noble-gas-functors-ms2/applicationLibrary/molecularDynamics/molecularDynamicsLibrary/KryptonPairFunctor.h
template <typename Particle_T>
class KryptonFunctorReference : public Functor<Particle_T> {
 public:
  explicit KryptonFunctorReference(
      double A, double a1, double a2, double a_m1, double b,
      double C6, double C8, double C10, double C12, double C14, double C16,
      bool newton3)
      : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b),
        _C6(C6), _C8(C8), _C10(C10),
        _C12(C12), _C14(C14), _C16(C16),
        _newton3(newton3) {}

  
  void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
    const auto& r1 = p1.getR();
    const auto& r2 = p2.getR();

    const double dx = r1[0] - r2[0];
    const double dy = r1[1] - r2[1];
    const double dz = r1[2] - r2[2];

    double r2sq = dx*dx + dy*dy + dz*dz;
    if (r2sq < 1e-24) r2sq = 1e-24;

    const double r    = std::sqrt(r2sq);


    const double expAlphaTerm = std::exp(_a1 * r + _a2 * r * r + _a_m1 / r);
    const double alphaTerm    = _a1 + 2.0 * _a2 * r - _a_m1 / (r * r);
    const double firstTerm    = (-_A * alphaTerm * expAlphaTerm) / r;

    const double bR = _b * r;
    const double expbr = std::exp(-bR);

    const auto TT = [&](int n, double C2n) {
      double inner = 0.0;
      for (int k = 0; k <= 2 * n; ++k)
        inner += std::pow(bR, k) / std::tgamma(k + 1.0); 
      double term = C2n * std::pow(r, -2.0 * n) *
                    (-2.0 * n + expbr * ((2.0 * n) * inner + bR * std::pow(bR, 2*n) / std::tgamma(2*n + 1.0)));
      return term;
    };

    const double term6  = TT(3, _C6);
    const double term8  = TT(4, _C8);
    const double term10 = TT(5, _C10);
    const double term12 = TT(6, _C12);
    const double term14 = TT(7, _C14);
    const double term16 = TT(8, _C16);

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16) / (r * r);
    const double Fmag = firstTerm + secondTerm;

    const double invr = 1.0 / r;
    const double coeff = Fmag * invr;

    std::array<double, 3> F{coeff * dx, coeff * dy, coeff * dz};


    p1.addF(F);
    if (_newton3) 
    p2.subF(F);

     
  }
    bool allowsNewton3() const { return true; }
    bool usesNewton3()   const  { return _newton3; }
 private:
  double _A, _a1, _a2, _a_m1, _b;
  double _C6, _C8, _C10, _C12, _C14, _C16;
  bool _newton3;
};
