#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "../Functors/Functor.h"
#include "../Particle.h"
#include "../generated_files/generated_krypton_force.hpp"
// Credits: https://github.com/AutoPas/AutoPas/blob/feat/3xa/noble-gas-functors-ms2/applicationLibrary/molecularDynamics/molecularDynamicsLibrary/KryptonPairFunctor.h
template <typename Particle_T>
class KryptonFunctorGenerated : public Functor<Particle_T> {
public:
  // Parameters: A, a1, a2, a_{-1}, b and C6, C8, C10
  // Higher dispersion coefficients C12, C14, C16 are derived internally
  // in krypton::computeForce via the recursion (Eq. 9).
  explicit KryptonFunctorGenerated(double A, double a1, double a2, double a_m1, double b,
                                   double C6, double C8, double C10,
                                   bool newton3 = true, double cutoff = 0.0)
      : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b),
        _C6(C6), _C8(C8), _C10(C10),
        _newton3(newton3), _cutoff(cutoff){}

  void AoSFunctor(Particle_T &p1, Particle_T &p2) override {
  const auto &r1 = p1.getR();
  const auto &r2 = p2.getR();

  const double dx = r1[0] - r2[0];
  const double dy = r1[1] - r2[1];
  const double dz = r1[2] - r2[2];

  constexpr double EPS = 1e-24;
  double r2sq = dx * dx + dy * dy + dz * dz;
  if (r2sq < EPS) r2sq = EPS;

  

  if (_cutoff > 0.0) {
    const double cutoff2 = _cutoff * _cutoff;
    if (r2sq > cutoff2) return;
  }

  const double r = std::sqrt(r2sq);

  double Fmag = krypton::computeForce(
      r, _A, _a1, _a2, _a_m1, _b, _C6, _C8, _C10);

  const double invr = 1.0 / r;
  const double coeff = Fmag * invr;

  std::array<double, 3> F{coeff * dx, coeff * dy, coeff * dz};

  p1.addF(F);
  if (_newton3) p2.subF(F);
}


  bool allowsNewton3() const  { return true; }
  bool usesNewton3() const  { return _newton3; }

private:
  double _A, _a1, _a2, _a_m1, _b;
  double _C6, _C8, _C10;
  bool _newton3;
  double _cutoff;
};
