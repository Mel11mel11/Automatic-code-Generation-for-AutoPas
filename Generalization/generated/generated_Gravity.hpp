
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <class Particle_T>
class GravityFunctor_Gen : public Functor<Particle_T> {
public:
    explicit GravityFunctor_Gen(bool newton3 = true, double G, double m1, double m2)
        : _newton3(newton3), _G(G), _m1(m1), _m2(m2) {}

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
        const double G = _G;
        const double m1 = _m1;
        const double m2 = _m2;

        // --- codegen: force magnitude F(r, params) ---
        // Fmag = -G*m1*m2/fast_pow(r, 2)
        const double Fmag = -G*m1*m2/fast_pow(r, 2);

        // Vector force: F = Fmag * r̂
        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double,3> F{fx, fy, fz};

        a.addF(F);
        if (_newton3) {
            b.subF(F);
        }
    }

private:
    bool _newton3;
    double _G;
    double _m1;
    double _m2;
};
