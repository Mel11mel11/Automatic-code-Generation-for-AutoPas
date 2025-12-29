#pragma once
#include "../ArrayMath.h"
#include "../ArrayUtils.h"
#include "Functor.h"
#include <cmath>

template <typename Particle_T>
class MieFunctorReference : public Functor<Particle_T> {
public:
    MieFunctorReference(
        double sigma,
        double epsilon,
        int n,
        int m,
        bool newton3 = true,
        double cutoff = 0.0
    )
        : _sigma(sigma),
          _epsilon(epsilon),
          _n(n),
          _m(m),
          _newton3(newton3),
          _cutoff(cutoff)
    {
        const double nd = static_cast<double>(_n);
        const double md = static_cast<double>(_m);
        // requires n > m
        _C = (nd / (nd - md)) * std::pow(nd / md, md / (nd - md));
    }

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {
        using namespace arrayMath::literals;

        const auto& r1 = p1.getR();
        const auto& r2 = p2.getR();
        auto dr = r1 - r2;

        double r2sq = arrayMath::dot(dr, dr);
        constexpr double minR2 = 1e-24;
        if (r2sq < minR2) r2sq = minR2;

        // ---- cutoff (AutoPas-style, local) ----
        if (_cutoff > 0.0) {
            const double cutoff2 = _cutoff * _cutoff;
            if (r2sq > cutoff2) return;
        }

        const double invr2 = 1.0 / r2sq;
        const double invr  = std::sqrt(invr2);

        const double sr   = _sigma * invr;
        const double sr_n = std::pow(sr, _n);
        const double sr_m = std::pow(sr, _m);

        // coefficient multiplying dr:
        // C*eps*( n*(σ/r)^n - m*(σ/r)^m ) / r^2
        const double coef =
            _C * _epsilon * (_n * sr_n - _m * sr_m) * invr2;

        auto F = dr * coef;
        p1.addF(F);
        if (_newton3) p2.subF(F);
    }

    bool usesNewton3() const  { return _newton3; }
    bool allowsNewton3() const  { return true; }

private:
    double _sigma;
    double _epsilon;
    double _C{};
    int    _n;
    int    _m;
    bool   _newton3;
    double _cutoff;
};
