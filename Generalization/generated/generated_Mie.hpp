
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>

template <class Particle_T>
class MieFunctor_Gen : public Functor<Particle_T> {
public:
    explicit MieFunctor_Gen(bool newton3 = true, double sigma, double epsilon, double n, double m, double C)
        : _newton3(newton3), _sigma(sigma), _epsilon(epsilon), _n(n), _m(m), _C(C) {}

    bool allowsNewton3() const override { return true; }
    bool usesNewton3()   const override { return _newton3; }

    void AoSFunctor(Particle_T& a, Particle_T& b) override {
        // Displacement a -> b (keep the same direction convention as reference)
        const auto& ra = a.getR();
        const auto& rb = b.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        // r^2 and r; guard against r -> 0
        constexpr double EPS = 1e-24;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;
        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        // Parameter aliases
        const double sigma = _sigma;
        const double epsilon = _epsilon;
        const double n = _n;
        const double m = _m;
        const double C = _C;

        // --- codegen: force magnitude F(r, params) ---
        // Fmag = -C*epsilon*(m*std::pow(sigma/r, m)/r - n*std::pow(sigma/r, n)/r)
        const double Fmag = -C*epsilon*(m*std::pow(sigma/r, m)/r - n*std::pow(sigma/r, n)/r);

        // Vector force: F = Fmag * r̂
        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        a.addF(fx, fy, fz);
        if (_newton3) {
            b.addF(-fx, -fy, -fz);
        }
    }

private:
    bool _newton3;
    double _sigma;
    double _epsilon;
    double _n;
    double _m;
    double _C;
};
