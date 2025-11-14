#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <class Particle_T>
class MieFunctor_Gen : public Functor<Particle_T> {
public:
    explicit MieFunctor_Gen(bool newton3,
                            double sigma,
                            double epsilon,
                            double n,
                            double m)
        : _newton3(newton3),
          _sigma(sigma),
          _epsilon(epsilon),
          _n(n),
          _m(m)
    {
        // Mie normalization constant (same as reference class)
        _C = (_n / (_n - _m)) * std::pow(_n / _m, _m / (_n - _m));
    }

    bool allowsNewton3() const  { return true; }
    bool usesNewton3()   const  { return _newton3; }

    void AoSFunctor(Particle_T& a, Particle_T& b) override {

        const auto& ra = a.getR();
        const auto& rb = b.getR();

        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        constexpr double EPS = 1e-24;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;

        const double r     = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        const double sr   = _sigma / r;
        const double sr_n = std::pow(sr, _n);
        const double sr_m = std::pow(sr, _m);

        // FORCE: F(r) = C ε ( n(σ/r)^n - m(σ/r)^m ) / r
        const double Fmag =
            _C * _epsilon * (_n * sr_n - _m * sr_m) / r;

        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double,3> F{fx, fy, fz};
        a.addF(F);
        if (_newton3) b.subF(F);
    }

private:
    bool   _newton3;

    double _sigma;
    double _epsilon;
    double _n;
    double _m;

    double _C; // now computed inside constructor
};
