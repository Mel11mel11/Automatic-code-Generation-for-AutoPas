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

template <typename Particle_T>
class KryptonFunctorReference : public Functor<Particle_T> {
 public:
  explicit KryptonFunctorReference(double A, double a1, double a2, double a_m1, double b,
                                   double C6, double C8, double C10,
                                   bool newton3)
      : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b),
        _C6(C6), _C8(C8), _C10(C10),
        _newton3(newton3) {}

  void AoSFunctor(Particle_T &p1, Particle_T &p2) override {
    const auto &r1 = p1.getR();
    const auto &r2 = p2.getR();

    const double dx = r1[0] - r2[0];
    const double dy = r1[1] - r2[1];
    const double dz = r1[2] - r2[2];

    double r2sq = dx * dx + dy * dy + dz * dz;
    if (r2sq < 1e-24) r2sq = 1e-24;

    const double r = std::sqrt(r2sq);

    // ----------------------------------------------------------
    // 1) Exponential repulsion term: -d/dr [ A exp(a1 r + a2 r² + a-1 / r) ]
    // ----------------------------------------------------------
    const double fexp = _a1 * r + _a2 * r * r + _a_m1 / r;
    const double expTerm = std::exp(fexp);

    const double dfdr = _a1 + 2.0 * _a2 * r - _a_m1 / (r * r);
    const double dV1dr = _A * dfdr * expTerm;

    // Force contribution (scalar radial magnitude)
    const double F1 = -dV1dr;

    // ----------------------------------------------------------
    // 2) Tang–Toennies damped dispersion (CORRECTED DERIVATIVE)
    //     V_n = - C2n / r^(2n) [1 - e^{-b r} S_{0..2n}(b r)]
    //     F_n = -dV_n/dr  --> (CORRECT FORM BELOW)
    // ----------------------------------------------------------

    const double bR = _b * r;
    const double expbr = std::exp(-bR);

    auto TT = [&](int n, double C2n) {
      //
      // Compute the polynomial sum S = sum_{k=0}^{2n} (bR)^k / k!
      //
      double S = 0.0;
      for (int k = 0; k <= 2 * n; ++k) {
        S += std::pow(bR, k) / std::tgamma(k + 1.0);
      }

      // Correct derivative of the Tang–Toennies-damped dispersion:
      //
      // F_n(r) =
      //     C2n * r^(-(2n+1)) *
      //       [ -2n + exp(-b r) * ( 2n * S + (bR)^(2n+1) / (2n)! ) ]
      //
      double term =
        C2n * std::pow(r, -(2 * n + 1.0)) *
        (
          -2.0 * n
          +
          expbr *
          (
            (2.0 * n) * S
            + std::pow(bR, 2 * n + 1) / std::tgamma(2.0 * n + 1.0)
          )
        );

      return term;   // already the correct radial force contribution
    };

    // Recursion for higher dispersion coefficients
    const double C12 = _C6  * std::pow(_C10 / _C8, 3.0);
    const double C14 = _C8  * std::pow(C12   / _C10, 3.0);
    const double C16 = _C10 * std::pow(C14   / C12, 3.0);

    // Force contributions from n=3..8 (C6..C16)
    const double F6  = TT(3, _C6);
    const double F8  = TT(4, _C8);
    const double F10 = TT(5, _C10);
    const double F12 = TT(6, C12);
    const double F14 = TT(7, C14);
    const double F16 = TT(8, C16);

    const double F2 = F6 + F8 + F10 + F12 + F14 + F16;

    // ----------------------------------------------------------
    // Total radial force
    // ----------------------------------------------------------
    const double Fmag = F1 + F2;

    const double invr = 1.0 / r;
    const double coeff = Fmag * invr;

    std::array<double, 3> F{coeff * dx, coeff * dy, coeff * dz};

    p1.addF(F);
    if (_newton3) {
      p2.subF(F);
    }
  }

  bool allowsNewton3() const  { return true; }
  bool usesNewton3() const    { return _newton3; }

 private:
  double _A, _a1, _a2, _a_m1, _b;
  double _C6, _C8, _C10;
  bool _newton3;
};
