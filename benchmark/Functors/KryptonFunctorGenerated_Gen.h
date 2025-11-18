#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <class Particle_T>
class KryptonFunctorGenerated_Gen : public Functor<Particle_T> {
public:
    explicit KryptonFunctorGenerated_Gen(double A, double a1, double a2, double a_m1,
                                         double b, double C6, double C8, double C10,
                                         bool newton3 = true)
        : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b),
          _C6(C6), _C8(C8), _C10(C10),
          _C12(0), _C14(0), _C16(0),
          _newton3(newton3)
    {
        // dispersion runtime compute
        _C12 = _C6 * std::pow(_C10 / _C8, 3.0);
        _C14 = _C8 * std::pow(_C12 / _C10, 3.0);
        _C16 = _C10 * std::pow(_C14 / _C12, 3.0);
    }

    void AoSFunctor(Particle_T& a, Particle_T& bPart) override {
        const auto& ra = a.getR();
        const auto& rb = bPart.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        constexpr double EPS = 1e-24;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;

        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        // Parameter aliases
        const double A  = _A;
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

        // --- YOUR GIANT FMAG FORMULA HERE ---
       const double Fmag = (1.0/20922789888000.0)*(
-20922789888000*A*fast_pow(C8, 15)*a1*fast_pow(r, 17)*std::exp(a1*r + a2*fast_pow(r, 2) + a_m1/r + b*r)
- 41845579776000*A*fast_pow(C8, 15)*a2*fast_pow(r, 18)*std::exp(a1*r + a2*fast_pow(r, 2) + a_m1/r + b*r)
+ 20922789888000*A*fast_pow(C8, 15)*a_m1*fast_pow(r, 15)*std::exp(a1*r + a2*fast_pow(r, 2) + a_m1/r + b*r)
+ fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 17)*fast_pow(r, 17)
+ 16*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 16)*fast_pow(r, 16)
+ 256*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 15)*fast_pow(r, 15)
+ 3840*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 14)*fast_pow(r, 14)
+ 53760*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 13)*fast_pow(r, 13)
+ 698880*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 12)*fast_pow(r, 12)
+ 8386560*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 11)*fast_pow(r, 11)
+ 92252160*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 10)*fast_pow(r, 10)
+ 922521600*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 9)*fast_pow(r, 9)
+ 8302694400*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 8)*fast_pow(r, 8)
+ 66421555200*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 7)*fast_pow(r, 7)
+ 464950886400*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 6)*fast_pow(r, 6)
+ 2789705318400*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 5)*fast_pow(r, 5)
+ 13948526592000*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 4)*fast_pow(r, 4)
+ 55794106368000*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 3)*fast_pow(r, 3)
+ 167382319104000*fast_pow(C10, 10)*fast_pow(C6, 6)*fast_pow(b, 2)*fast_pow(r, 2)
+ 334764638208000*fast_pow(C10, 10)*fast_pow(C6, 6)*b*r
- 334764638208000*fast_pow(C10, 10)*fast_pow(C6, 6)*std::exp(b*r)
+ 334764638208000*fast_pow(C10, 10)*fast_pow(C6, 6)
+ 240*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 15)*fast_pow(r, 17)
+ 3360*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 14)*fast_pow(r, 16)
+ 47040*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 13)*fast_pow(r, 15)
+ 611520*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 12)*fast_pow(r, 14)
+ 7338240*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 11)*fast_pow(r, 13)
+ 80720640*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 10)*fast_pow(r, 12)
+ 807206400*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 9)*fast_pow(r, 11)
+ 7264857600*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 8)*fast_pow(r, 10)
+ 58118860800*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 7)*fast_pow(r, 9)
+ 406832025600*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 6)*fast_pow(r, 8)
+ 2440992153600*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 5)*fast_pow(r, 7)
+ 12204960768000*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 4)*fast_pow(r, 6)
+ 48819843072000*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 3)*fast_pow(r, 5)
+ 146459529216000*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(b, 2)*fast_pow(r, 4)
+ 292919058432000*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*b*fast_pow(r, 3)
- 292919058432000*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(r, 2)*std::exp(b*r)
+ 292919058432000*fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7)*fast_pow(r, 2)
+ 43680*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 13)*fast_pow(r, 17)
+ 524160*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 12)*fast_pow(r, 16)
+ 6289920*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 11)*fast_pow(r, 15)
+ 69189120*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 10)*fast_pow(r, 14)
+ 691891200*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 9)*fast_pow(r, 13)
+ 6227020800*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 8)*fast_pow(r, 12)
+ 49816166400*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 7)*fast_pow(r, 11)
+ 348713164800*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 6)*fast_pow(r, 10)
+ 2092278988800*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 5)*fast_pow(r, 9)
+ 10461394944000*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 4)*fast_pow(r, 8)
+ 41845579776000*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 3)*fast_pow(r, 7)
+ 125536739328000*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(b, 2)*fast_pow(r, 6)
+ 251073478656000*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*b*fast_pow(r, 5)
- 251073478656000*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(r, 4)*std::exp(b*r)
+ 251073478656000*fast_pow(C10, 3)*C6*fast_pow(C8, 12)*fast_pow(r, 4)
+ 5765760*C10*fast_pow(C8, 15)*fast_pow(b, 11)*fast_pow(r, 17)
+ 57657600*C10*fast_pow(C8, 15)*fast_pow(b, 10)*fast_pow(r, 16)
+ 576576000*C10*fast_pow(C8, 15)*fast_pow(b, 9)*fast_pow(r, 15)
+ 5189184000*C10*fast_pow(C8, 15)*fast_pow(b, 8)*fast_pow(r, 14)
+ 41513472000*C10*fast_pow(C8, 15)*fast_pow(b, 7)*fast_pow(r, 13)
+ 290594304000*C10*fast_pow(C8, 15)*fast_pow(b, 6)*fast_pow(r, 12)
+ 1743565824000*C10*fast_pow(C8, 15)*fast_pow(b, 5)*fast_pow(r, 11)
+ 8717829120000*C10*fast_pow(C8, 15)*fast_pow(b, 4)*fast_pow(r, 10)
+ 34871316480000*C10*fast_pow(C8, 15)*fast_pow(b, 3)*fast_pow(r, 9)
+ 104613949440000*C10*fast_pow(C8, 15)
+ 209227898880000*C10*fast_pow(C8, 15)*b*fast_pow(r, 7)
- 209227898880000*C10*fast_pow(C8, 15)*fast_pow(r, 6)*std::exp(b*r)
+ 209227898880000*C10*fast_pow(C8, 15)*fast_pow(r, 6)
+ 29059430400*C6*fast_pow(C8, 15)*fast_pow(b, 7)*fast_pow(r, 17)
+ 174356582400*C6*fast_pow(C8, 15)*fast_pow(b, 6)*fast_pow(r, 16)
+ 1046139494400*C6*fast_pow(C8, 15)*fast_pow(b, 5)*fast_pow(r, 15)
+ 5230697472000*C6*fast_pow(C8, 15)*fast_pow(b, 4)*fast_pow(r, 14)
+ 20922789888000*C6*fast_pow(C8, 15)*fast_pow(b, 3)*fast_pow(r, 13)
+ 62768369664000*C6*fast_pow(C8, 15)*fast_pow(b, 2)*fast_pow(r, 12)
+ 125536739328000*C6*fast_pow(C8, 15)*b*fast_pow(r, 11)
- 125536739328000*C6*fast_pow(C8, 15)*fast_pow(r, 10)*std::exp(b*r)
+ 125536739328000*C6*fast_pow(C8, 15)*fast_pow(r, 10)
+ 518918400*fast_pow(C8, 16)*fast_pow(b, 9)*fast_pow(r, 17)
+ 4151347200*fast_pow(C8, 16)*fast_pow(b, 8)*fast_pow(r, 16)
+ 33210777600*fast_pow(C8, 16)*fast_pow(b, 7)*fast_pow(r, 15)
+ 232475443200*fast_pow(C8, 16)*fast_pow(b, 6)*fast_pow(r, 14)
+ 1394852659200*fast_pow(C8, 16)*fast_pow(b, 5)*fast_pow(r, 13)
+ 6974263296000*fast_pow(C8, 16)*fast_pow(b, 4)*fast_pow(r, 12)
+ 27897053184000*fast_pow(C8, 16)*fast_pow(b, 3)*fast_pow(r, 11)
+ 83691159552000*fast_pow(C8, 16)*fast_pow(b, 2)*fast_pow(r, 10)
+ 167382319104000*fast_pow(C8, 16)*b*fast_pow(r, 9)
- 167382319104000*fast_pow(C8, 16)*fast_pow(r, 8)*std::exp(b*r)
+ 167382319104000*fast_pow(C8, 16)*fast_pow(r, 8)
)*std::exp(-b*r)/(fast_pow(C8, 15)*fast_pow(r, 17));


        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double,3> F{fx, fy, fz};
        a.addF(F);
        if (_newton3) {
            bPart.subF(F);
        }
    }

    bool allowsNewton3() const { return true; }
    bool usesNewton3()   const { return _newton3; }

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
};
