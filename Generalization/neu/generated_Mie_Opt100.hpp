
#pragma once


#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>

#ifdef USE_FAST_POW
#include "FastPow.hpp"
#endif


template <class Particle_T>
class MieFunctor_Gen_Opt100 : public Functor<Particle_T> {
public:
    explicit MieFunctor_Gen_Opt100(double sigma, double epsilon, double n, double m, bool newton3 = true, double cutoff = 0.0)
        : _sigma(sigma), _epsilon(epsilon), _n(n), _m(m), _C(0.0), _newton3(newton3), _cutoff(cutoff)
    {

        _C = (_n / (_n - _m)) * std::pow(_n / _m, _m / (_n - _m));


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
        const double n = _n;
        const double m = _m;
        const double C = _C;

        #ifdef USE_MASS
        const double p1m = p1.getMass();
        const double p2m = p2.getMass();
        #endif




        const double Fmag = -C*epsilon*inv_r*(m*std::pow(inv_r*sigma, m) - n*std::pow(inv_r*sigma, n));

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
    double _n;
    double _m;
    double _C;
    bool _newton3;
    double _cutoff;
};
