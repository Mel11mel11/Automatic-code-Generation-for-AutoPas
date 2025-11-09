#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"

struct KryptonParams {
  double alpha1;      // _alpha1
  double alpha2;      // _alpha2
  double alphaInv1;   // _alphaInv1
  double constA;      // _constA
  double constATilde; // _constATilde
  double alphaTilde;  // _alphaTilde
  double constb;      // _constb
  double C6;  double C8;  double C10;
  double C12; double C14; double C16;
  double minDistance; // _minDistance
};

template <typename Particle_T>
class KryptonFunctorReference : public Functor<Particle_T> {
 public:
  explicit KryptonFunctorReference(double cutoff, KryptonParams p)
      : _cutoffSquared{cutoff * cutoff},
        _p{std::move(p)} {}

  bool allowsNewton3() const override { return true; }

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {
    using namespace autopas::utils::ArrayMath::literals;

    const auto displacement = i.getR() - j.getR();
    const double dist2 = autopas::utils::ArrayMath::dot(displacement, displacement);
    if (dist2 > _cutoffSquared) return;

    // Precalculations
    const double dist    = std::sqrt(dist2);
    const double distInv = 1.0 / (dist + 1e-300);  // r->0 smoothing 
    const double distInv2 = distInv * distInv;

    const double distInv6  = distInv2 * distInv2 * distInv2;  // r^-6
    const double distNeg8  = distInv6 * distInv2;             // r^-8
    const double distNeg10 = distNeg8 * distInv2;             // r^-10
    const double distNeg12 = distNeg10 * distInv2;            // r^-12
    const double distNeg14 = distNeg12 * distInv2;            // r^-14
    const double distNeg16 = distNeg14 * distInv2;            // r^-16

    
    const double expAlphaTerm = std::exp(_p.alpha1 * dist + _p.alpha2 * dist2 + _p.alphaInv1 * distInv);
   
    const double alphaTerm = _p.alpha1 + 2.0 * _p.alpha2 * dist - _p.alphaInv1 * distInv2;

    
    const double firstTerm = (-_p.constA * alphaTerm * expAlphaTerm) * distInv;

    
    const double bdist = _p.constb * dist;
    
    const auto &F = _invFactorials;

    
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

    const double ksum0  = 1.0;
    const double ksum1  = bdist;
    const double ksum2  = b2  * F[2];
    const double ksum3  = b3  * F[3];
    const double ksum4  = b4  * F[4];
    const double ksum5  = b5  * F[5];
    const double ksum6  = b6  * F[6];
    const double ksum7  = b7  * F[7];
    const double ksum8  = b8  * F[8];
    const double ksum9  = b9  * F[9];
    const double ksum10 = b10 * F[10];
    const double ksum11 = b11 * F[11];
    const double ksum12 = b12 * F[12];
    const double ksum13 = b13 * F[13];
    const double ksum14 = b14 * F[14];
    const double ksum15 = b15 * F[15];
    const double ksum16 = b16 * F[16];

    const double ksumacc6  = ksum0 + ksum1 + ksum2 + ksum3 + ksum4 + ksum5 + ksum6;
    const double ksumacc8  = ksumacc6  + ksum7 + ksum8;
    const double ksumacc10 = ksumacc8  + ksum9 + ksum10;
    const double ksumacc12 = ksumacc10 + ksum11 + ksum12;
    const double ksumacc14 = ksumacc12 + ksum13 + ksum14;
    const double ksumacc16 = ksumacc14 + ksum15 + ksum16;

    const double expbr = std::exp(-bdist);

    // Cn expressions
    const double term6  = _p.C6  * distInv6  * (-6.0  + expbr * ( 6.0  * ksumacc6  + bdist * ksum6 ));
    const double term8  = _p.C8  * distNeg8  * (-8.0  + expbr * ( 8.0  * ksumacc8  + bdist * ksum8 ));
    const double term10 = _p.C10 * distNeg10 * (-10.0 + expbr * (10.0 * ksumacc10 + bdist * ksum10));
    const double term12 = _p.C12 * distNeg12 * (-12.0 + expbr * (12.0 * ksumacc12 + bdist * ksum12));
    const double term14 = _p.C14 * distNeg14 * (-14.0 + expbr * (14.0 * ksumacc14 + bdist * ksum14));
    const double term16 = _p.C16 * distNeg16 * (-16.0 + expbr * (16.0 * ksumacc16 + bdist * ksum16));

    const double secondTerm = (term6 + term8 + term10 + term12 + term14 + term16) * distInv2;

    
    const double scalar =
        (dist >= _p.minDistance)
            ? (firstTerm + secondTerm)
            : (_p.constATilde * distInv2 * std::exp(-_p.alphaTilde * dist) * (distInv + _p.alphaTilde));

    const auto f = displacement * scalar;

    i.addF(f);
    if (newton3) j.subF(f);
  }

 private:
  double _cutoffSquared;

  KryptonParams _p;

  static constexpr std::array<double, 17> _invFactorials = {
      1.0,
      1.0,
      1.0 / 2.0,
      1.0 / 6.0,
      1.0 / 24.0,
      1.0 / 120.0,
      1.0 / 720.0,
      1.0 / 5040.0,
      1.0 / 40320.0,
      1.0 / 362880.0,
      1.0 / 3628800.0,
      1.0 / 39916800.0,
      1.0 / 479001600.0,
      1.0 / 6227020800.0,
      1.0 / 87178291200.0,
      1.0 / 1307674368000.0,
      1.0 / 20922789888000.0};
};
template <typename Particle_T>
constexpr std::array<double,17> KryptonFunctorReference<Particle_T>::_invFactorials;