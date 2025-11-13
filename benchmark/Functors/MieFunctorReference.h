
#pragma once
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include <cmath>

template <typename Particle_T>
class MieFunctorReference : public Functor<Particle_T> {
public:
    MieFunctorReference(double sigma, double epsilon, int n, int m, bool newton3 = true)
        : _sigma(sigma), _epsilon(epsilon), _n(n), _m(m), _newton3(newton3) {
        const double nd = static_cast<double>(_n);
        const double md = static_cast<double>(_m);
        // requires n > m
        _C = (nd / (nd - md)) * std::pow(nd / md, md / (nd - md));
    }

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        using namespace arrayMath::literals;  // enables +,-,*, dot for std::array

        const auto& r1 = p1.getR();
        const auto& r2 = p2.getR();
        auto dr = r1 - r2;                         
        double r2sq = arrayMath::dot(dr, dr);
        if (r2sq < 1e-24) r2sq = 1e-24;            
        const double invr2 = 1.0 / r2sq;
        const double invr  = std::sqrt(invr2);     

        const double sr   = _sigma * invr;        
        const double sr_n = std::pow(sr, _n);
        const double sr_m = std::pow(sr, _m);

        // coefficient multiplying dr: C*eps*( n*(σ/r)^n - m*(σ/r)^m ) / r^2
        const double coef = _C * _epsilon * (_n * sr_n - _m * sr_m) * invr2;

        auto F = dr * coef;                         // vector force on p1
        p1.addF(F);
        if (_newton3) p2.subF(F);
    }
    bool usesNewton3() const  { return _newton3; }
    bool allowsNewton3() const { return true; }

private:
    double _sigma, _epsilon, _C{};
    int _n, _m;
    bool _newton3;
};
