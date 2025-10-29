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

        auto dr = p1.getR() - p2.getR();          // vector with 3 dimensions x,y,z
        double dr2 = arrayMath::dot(dr, dr);      // squared distance and now it is scalar
        constexpr double minR2 = 1e-8;     // üstte tanımla istersen
        if (dr2 < minR2) dr2 = minR2;      // <-- değişken adı: dr2
           // avoid division by 0 yeah yeah that is a must

        double invdr2 = 1.0 / dr2;                  // that is for 1/r^2
        double lj6 = _sigmaSquared * invdr2;       // (sigma^2 / r^2)
        lj6 = lj6 * lj6 * lj6;                    // (sigma^6 / r^6)
        double lj12 = lj6 * lj6;                  // (sigma^12 / r^12)
        double lj12m6 = lj12 - lj6;               // (sigma^12 - sigma^6)
        double fac = _epsilon24 * (lj12 + lj12m6) * invdr2; // dV/dr factor
        auto F = dr * fac;   // multiply vector by scalar to get force vector

        p1.addF(F);
        if (_newton3) p2.subF(F);                 
    }

private:
    double _sigmaSquared{};
    double _epsilon24{};
    bool _newton3;                                
};
