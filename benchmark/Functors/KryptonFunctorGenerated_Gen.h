
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>

template<class Particle_T>
class KryptonFunctorGenerated_Gen : public Functor<Particle_T> {
public:
    explicit KryptonFunctorGenerated_Gen(double A, double a1, double a2, double a_m1, double b, double C6, double C8, double C10, bool newton3 = true)
        : _A(A), _a1(a1), _a2(a2), _a_m1(a_m1), _b(b), _C6(C6), _C8(C8), _C10(C10), _C12(0.0), _C14(0.0), _C16(0.0), _newton3(newton3)
    {
        // Precompute higher dispersion coefficients once per functor
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
        if(r2 < EPS) r2 = EPS;

        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

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

        const double br = b * r;
        const double exp_br = std::exp(-br);

        // powers of 1/r
        const double inv_r2 = 1.0 / (r * r);
        const double inv_r4  = inv_r2 * inv_r2;
        const double inv_r6  = inv_r4 * inv_r2;
        const double inv_r8  = inv_r6 * inv_r2;
        const double inv_r10 = inv_r8 * inv_r2;
        const double inv_r12 = inv_r10 * inv_r2;
        const double inv_r14 = inv_r12 * inv_r2;
        const double inv_r16 = inv_r14 * inv_r2;
        const double inv_r18 = inv_r16 * inv_r2;

        // Unrolled series for Sum_{k=0}^n (br^k / k!)
        double t1  = br;
        double t2  = t1 * br / 2.0;
        double t3  = t2 * br / 3.0;
        double t4  = t3 * br / 4.0;
        double t5  = t4 * br / 5.0;
        double t6  = t5 * br / 6.0;
        double t7  = t6 * br / 7.0;
        double t8  = t7 * br / 8.0;
        double t9  = t8 * br / 9.0;
        double t10 = t9 * br / 10.0;
        double t11 = t10 * br / 11.0;
        double t12 = t11 * br / 12.0;
        double t13 = t12 * br / 13.0;
        double t14 = t13 * br / 14.0;
        double t15 = t14 * br / 15.0;
        double t16 = t15 * br / 16.0;

        double S1  = 1.0 + t1;
        double S2  = S1  + t2;
        double S3  = S2  + t3;
        double S4  = S3  + t4;
        double S5  = S4  + t5;
        double S6  = S5  + t6;
        double S7  = S6  + t7;
        double S8  = S7  + t8;
        double S9  = S8  + t9;
        double S10 = S9  + t10;
        double S11 = S10 + t11;
        double S12 = S11 + t12;
        double S13 = S12 + t13;
        double S14 = S13 + t14;
        double S15 = S14 + t15;
        double S16 = S15 + t16;

        const double D6  = 1.0 - exp_br * S6;
        const double D8  = 1.0 - exp_br * S8;
        const double D10 = 1.0 - exp_br * S10;
        const double D12 = 1.0 - exp_br * S12;
        const double D14 = 1.0 - exp_br * S14;
        const double D16 = 1.0 - exp_br * S16;

        const double exp_main = std::exp(a1*r + a2*r*r + a_m1/r + b*r);

        double F = A * exp_main * (a1 + 2.0*a2*r - a_m1/(r*r) + b);

        F += C6  * (-6.0  * inv_r8 ) * D6;
        F += C8  * (-8.0  * inv_r10) * D8;
        F += C10 * (-10.0 * inv_r12) * D10;
        F += C12 * (-12.0 * inv_r14) * D12;
        F += C14 * (-14.0 * inv_r16) * D14;
        F += C16 * (-16.0 * inv_r18) * D16;

        const double Fmag = F;

        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double,3> Fv{{fx, fy, fz}};
        p1.addF(Fv);
        if(_newton3) p2.subF(Fv);
    }

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
