#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include "../ArrayMath.h"   
#include <cmath>

template <class Particle_T>
class Automated_LJ : public Functor<Particle_T> {
public:
    
    Automated_LJ(double sigma, double epsilon, bool newton3 = true)
        : _newton3(newton3), _sigma(sigma), _epsilon(epsilon) {}

    bool allowsNewton3() const { return true; }
    bool usesNewton3()  const { return _newton3; }

    void AoSFunctor(Particle_T& a, Particle_T& b) override {
        const auto& ra = a.getR();
        const auto& rb = b.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        constexpr double EPS = 1e-24;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;
        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        const double sigma   = _sigma;
        const double epsilon = _epsilon;

        // F magnitude (skaler)
        const double Fmag =
            -4.0 * epsilon * ( 6.0 * std::pow(sigma, 6)  / std::pow(r, 7)
                             - 12.0 * std::pow(sigma, 12) / std::pow(r,13) );

        // Bileşenler
        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        // >>> fark burada: tek bir F vektörü oluştur ve onu ver
        arrayMath::Array3D<double> F = {fx, fy, fz};

        a.addF(F);
        if (_newton3) {
            b.subF(F);   // veya: arrayMath::Array3D<double> Fn = {-fx,-fy,-fz}; b.addF(Fn);
        }
    }

private:
    bool _newton3;
    double _sigma;
    double _epsilon;
};
