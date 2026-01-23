
#pragma once


#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>

#ifdef USE_FAST_POW
#include "FastPow.hpp"
#endif


template <class Particle_T>
class KryptonFunctorGenerated_Gen_Opt010 : public Functor<Particle_T> {
public:
    explicit KryptonFunctorGenerated_Gen_Opt010(double A, double a1, double a2, double a_m1, double b, double C6, double C8, double C10, bool newton3 = true, double cutoff = 0.0)
        : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b), _C6(C6), _C8(C8), _C10(C10), _C12(0.0), _C14(0.0), _C16(0.0), _newton3(newton3), _cutoff(cutoff)
    {


        _C12 = _C6 * std::pow(_C10 / _C8, 3.0);
        _C14 = _C8 * std::pow(_C12 / _C10, 3.0);
        _C16 = _C10 * std::pow(_C14 / _C12, 3.0);

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
        const double A = _A;
        const double a1 = _a1;
        const double a2 = _a2;
        const double a_m1 = _a_m1;
        const double b = _b;
        const double C6 = _C6;
        const double C8 = _C8;
        const double C10 = _C10;
        const double C12 = _C12;
        const double C14 = _C14;
        const double C16 = _C16;

        #ifdef USE_MASS
        const double p1m = p1.getMass();
        const double p2m = p2.getMass();
        #endif


        const double tt_x = b*r;
        double tt_arr[17];
        tt_arr[0] = 1.0;
        double tt_pow = 1.0;
        double tt_inv_fact = 1.0;
        double tt_sum = 1.0;
        for (int k = 1; k <= 16; ++k) {
          tt_pow *= tt_x;
          tt_inv_fact /= double(k);
          tt_sum += tt_pow * tt_inv_fact;
          tt_arr[k] = tt_sum;
        }
        const double tt5 = tt_arr[5];
        const double tt6 = tt_arr[6];
        const double tt7 = tt_arr[7];
        const double tt8 = tt_arr[8];
        const double tt9 = tt_arr[9];
        const double tt10 = tt_arr[10];
        const double tt11 = tt_arr[11];
        const double tt12 = tt_arr[12];
        const double tt13 = tt_arr[13];
        const double tt14 = tt_arr[14];
        const double tt15 = tt_arr[15];
        const double tt16 = tt_arr[16];
        const double x0 = std::exp(-b*r);
        const double x1 = b*x0;
        const double x2 = C6*x1;
        const double x3 = std::pow(C10, 3)/std::pow(C8, 3);
        const double x4 = std::pow(C10, 6)*std::pow(C6, 3)/std::pow(C8, 8);
        const double x5 = std::pow(C10, 10)*std::pow(C6, 6)/std::pow(C8, 15);
        const double x6 = std::pow(r, 2);

        const double Fmag = -A*(a1 + 2*a2*r - a_m1/x6)*std::exp(a1*r + a2*x6 + a_m1/r) + C10*x1*(tt10 - tt9)/std::pow(r, 10) + 10*C10*(tt10*x0 - 1)/std::pow(r, 11) + 6*C6*(tt6*x0 - 1)/std::pow(r, 7) + 12*C6*x3*(tt12*x0 - 1)/std::pow(r, 13) - C8*x1*(tt7 - tt8)/std::pow(r, 8) + 8*C8*(tt8*x0 - 1)/std::pow(r, 9) - x2*(tt5 - tt6)/std::pow(r, 6) - x2*x3*(tt11 - tt12)/std::pow(r, 12) - x1*x4*(tt13 - tt14)/std::pow(r, 14) + 14*x4*(tt14*x0 - 1)/std::pow(r, 15) - x1*x5*(tt15 - tt16)/std::pow(r, 16) + 16*x5*(tt16*x0 - 1)/std::pow(r, 17);

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
    double _A;
    double _a1;
    double _a2;
    double _a_m1;
    double _b;
    double _C6;
    double _C8;
    double _C10;
    double _C12;
    double _C14;
    double _C16;
    bool _newton3;
    double _cutoff;
};
