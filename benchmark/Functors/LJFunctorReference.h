#pragma once
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include <iostream>

template <typename Particle_T>
class LJFunctorReference : public Functor<Particle_T> {
public:
    // Constructor with optional Newton3 flag
    LJFunctorReference(double sigma, double epsilon, bool newton3 = false)
        : Functor<Particle_T>(),
          _sigmaSquared(sigma * sigma),
          _epsilon24(epsilon * 24.0),
          _newton3(newton3) {}  // ✅ burayı ekledik

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        using namespace arrayMath::literals;

        auto dr = p1.getR() - p2.getR();          // displacement vector
        double dr2 = arrayMath::dot(dr, dr);      // squared distance
        if (dr2 < 1e-24) dr2 = 1e-24;             // avoid division by 0

        double invdr2 = 1.0 / dr2;
        double lj6 = _sigmaSquared * invdr2;
        lj6 = lj6 * lj6 * lj6;                    // (sigma^6 / r^6)
        double lj12 = lj6 * lj6;                  // (sigma^12 / r^12)
        double lj12m6 = lj12 - lj6;               // (sigma^12 - sigma^6)
        double fac = _epsilon24 * (lj12 + lj12m6) * invdr2; // dV/dr factor

        auto f = dr * fac;                       
        p1.addF(f);
        if (_newton3) p2.subF(f);                 
    }

private:
    double _sigmaSquared{};
    double _epsilon24{};
    bool _newton3;                                
};
