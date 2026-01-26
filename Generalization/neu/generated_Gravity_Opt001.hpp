
#pragma once
#define USE_FAST_POW

#define USE_MASS

#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>

#ifdef USE_FAST_POW
#include "FastPow.hpp"
#endif


template <class Particle_T>
class GravityFunctor_Gen_Opt001 : public Functor<Particle_T> {
public:
    explicit GravityFunctor_Gen_Opt001(double G, bool newton3 = true, double cutoff = 0.0)
        : _G(G), _newton3(newton3), _cutoff(cutoff)
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
        const double G = _G;

        #ifdef USE_MASS
        const double p1m = p1.getMass();
        const double p2m = p2.getMass();
        #endif




        const double Fmag = -G*fast_pow(inv_r, 2)*p1m*p2m;

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
    double _G;
    bool _newton3;
    double _cutoff;
};
