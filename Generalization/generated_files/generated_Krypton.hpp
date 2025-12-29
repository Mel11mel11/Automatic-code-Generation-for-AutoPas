
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <class Particle_T>
class KryptonFunctorGenerated_Gen : public Functor<Particle_T> {
public:
    explicit KryptonFunctorGenerated_Gen(double A, double a1, double a2, double a_m1, double b, double C6, double C8, double C10, bool newton3 = true, double cutoff = 0.0)
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

        const double p1m = p1.getMass();
        const double p2m = p2.getMass();

        const double x0 = fast_pow(inv_r, 11);
        const double x1 = fast_pow(inv_r, 7);
        const double x2 = fast_pow(inv_r, 9);
        const double x3 = fast_pow(inv_r, 13);
        const double x4 = fast_pow(C10, 3);
        const double x5 = fast_pow(C8, -3);
        const double x6 = fast_pow(inv_r, 15);
        const double x7 = fast_pow(C10, 6);
        const double x8 = fast_pow(C6, 3);
        const double x9 = fast_pow(C8, -8);
        const double x10 = fast_pow(C10, 10);
        const double x11 = fast_pow(C6, 6);
        const double x12 = fast_pow(C8, -15);
        const double x13 = fast_pow(inv_r, -1);
        const double x14 = fast_pow(inv_r, 2);
        const double x15 = fast_pow(x14, -1);
        const double x16 = std::exp(a1*x13 + a2*x15 + a_m1*inv_r);
        const double x17 = A*x16;
        const double x18 = fast_pow(inv_r, 6);
        const double x19 = b*x13;
        const double x20 = std::exp(-x19);
        const double x21 = 720*b;
        const double x22 = fast_pow(b, 7);
        const double x23 = fast_pow(x18, -1);
        const double x24 = b + 5*inv_r;
        const double x25 = b + inv_r;
        const double x26 = fast_pow(b, 5)/fast_pow(inv_r, 5);
        const double x27 = b + 4*inv_r;
        const double x28 = fast_pow(b, 2)*x15;
        const double x29 = b + 3*inv_r;
        const double x30 = fast_pow(b, 3)/fast_pow(inv_r, 3);
        const double x31 = b + 2*inv_r;
        const double x32 = fast_pow(b, 4)/fast_pow(inv_r, 4);
        const double x33 = fast_pow(inv_r, 8);
        const double x34 = fast_pow(b, 9);
        const double x35 = fast_pow(x33, -1);
        const double x36 = b + 7*inv_r;
        const double x37 = x22/x1;
        const double x38 = b + 6*inv_r;
        const double x39 = fast_pow(b, 6)*x23;
        const double x40 = fast_pow(inv_r, 10);
        const double x41 = fast_pow(b, 11);
        const double x42 = fast_pow(x40, -1);
        const double x43 = b + 9*inv_r;
        const double x44 = x34/x2;
        const double x45 = b + 8*inv_r;
        const double x46 = fast_pow(b, 8)*x35;
        const double x47 = fast_pow(inv_r, 12);
        const double x48 = fast_pow(b, 13);
        const double x49 = fast_pow(x47, -1);
        const double x50 = b + 11*inv_r;
        const double x51 = x41/x0;
        const double x52 = b + 10*inv_r;
        const double x53 = fast_pow(b, 10)*x42;
        const double x54 = fast_pow(inv_r, 14);
        const double x55 = fast_pow(b, 15);
        const double x56 = fast_pow(x54, -1);
        const double x57 = b + 13*inv_r;
        const double x58 = x48/x3;
        const double x59 = b + 12*inv_r;
        const double x60 = fast_pow(b, 12)*x49;
        const double x61 = fast_pow(inv_r, 16);

        const double Fmag = A*a_m1*x14*x16 - 10*C10*x0 + (1.0/3628800.0)*C10*x20*x40*(3628800*b + 36288000*inv_r + 3628800*x19*x43 + 30240*x24*x26 + 10*x25*x44 + 5040*x27*x39 + 1814400*x28*x45 + 720*x29*x37 + 604800*x30*x36 + 90*x31*x46 + 151200*x32*x38 + x41*x42) - 6*C6*x1 + (1.0/720.0)*C6*x18*x20*(4320*inv_r + x13*x21*x24 + x21 + x22*x23 + 6*x25*x26 + 360*x27*x28 + 120*x29*x30 + 30*x31*x32) + (1.0/479001600.0)*C6*x20*x4*x47*x5*(479001600*b + 5748019200*inv_r + 479001600*x19*x50 + 95040*x24*x37 + 12*x25*x51 + 3991680*x26*x36 + 11880*x27*x46 + 239500800*x28*x52 + 1320*x29*x44 + 79833600*x30*x43 + 132*x31*x53 + 19958400*x32*x45 + 665280*x38*x39 + x48*x49) - 12*C6*x3*x4*x5 - 8*C8*x2 + (1.0/40320.0)*C8*x20*x33*(40320*b + 322560*inv_r + 40320*x19*x36 + 6720*x24*x30 + 8*x25*x37 + 336*x26*x29 + 1680*x27*x32 + 20160*x28*x38 + 56*x31*x39 + x34*x35) - a1*x17 - 2*a2*x13*x17 - 16*fast_pow(inv_r, 17)*x10*x11*x12 + (1.0/20922789888000.0)*x10*x11*x12*x20*x61*(fast_pow(b, 17)/x61 + 240*fast_pow(b, 14)*x31*x56 + 20922789888000*b + 334764638208000*inv_r + 20922789888000*x19*(b + 15*inv_r) + 524160*x24*x51 + 16*x25*x55/x6 + 174356582400*x26*x50 + 43680*x27*x60 + 10461394944000*x28*(b + 14*inv_r) + 3360*x29*x58 + 3487131648000*x30*x57 + 871782912000*x32*x59 + 57657600*x36*x44 + 4151347200*x37*x43 + 5765760*x38*x53 + 29059430400*x39*x52 + 518918400*x45*x46) + (1.0/87178291200.0)*x20*x54*x7*x8*x9*(87178291200*b + 1220496076800*inv_r + 87178291200*x19*x57 + 240240*x24*x44 + 14*x25*x58 + 726485760*x26*x43 + 24024*x27*x53 + 43589145600*x28*x59 + 2184*x29*x51 + 14529715200*x30*x50 + 182*x31*x60 + 3632428800*x32*x52 + 17297280*x36*x37 + 2162160*x38*x46 + 121080960*x39*x45 + x55*x56) - 14*x6*x7*x8*x9;

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
