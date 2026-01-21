
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <class Particle_T>
class LJFunctor_Gen_O011 : public Functor<Particle_T> {
public:
    explicit LJFunctor_Gen_O011(double sigma, double epsilon, bool newton3 = true, double cutoff = 0.0)
        : _sigma(sigma), _epsilon(epsilon), _newton3(newton3), _cutoff(cutoff)
    {


    }

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        const auto& ra = p1.getR();
        const auto& rb = p2.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        constexpr double EPS = 1e-24;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;
        const double cutoff = _cutoff;
        const double cutoff2 = cutoff * cutoff;
        if (cutoff > 0.0 && r2 > cutoff2) return;
        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        // Parameter aliases
        const double sigma = _sigma;
        const double epsilon = _epsilon;

      
        const double x0 = fast_pow(sigma, 6);

        const double Fmag = 24*epsilon*fast_pow(inv_r, 7)*x0*(2*fast_pow(inv_r, 6)*x0 - 1);

        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double,3> F{fx, fy, fz};
        p1.addF(F);
        if (_newton3) p2.subF(F);
    }

    bool allowsNewton3() const { return true; }
    bool usesNewton3() const { return _newton3; }

private:
    double _sigma;
    double _epsilon;
    bool _newton3;
    double _cutoff;
};
