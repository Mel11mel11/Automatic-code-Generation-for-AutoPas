#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>

template <class SoAView_T>
class KryptonFunctorGenerated_Gen_SoA {
public:
    using SoAView = SoAView_T;

    explicit KryptonFunctorGenerated_Gen_SoA(
        double A,
        double a1,
        double a2,
        double a_m1,
        double b,
        double C6,
        double C8,
        double C10,
        bool newton3 = true
    )
      : _newton3(newton3)
      , _A(A)
      , _a1(a1)
      , _a2(a2)
      , _a_m1(a_m1)
      , _b(b)
      , _C6(C6)
      , _C8(C8)
      , _C10(C10)
    {
    }

    // =============================
    //        SOA FUNCTOR
    // =============================
    void SoAFunctor(SoAView &soa) {

        double *x  = soa.x.data();
        double *y  = soa.y.data();
        double *z  = soa.z.data();

        double *fx = soa.fx.data();
        double *fy = soa.fy.data();
        double *fz = soa.fz.data();

        double *mass = soa.mass.data();

        const std::size_t N = soa.size();

        for (std::size_t i = 0; i < N; ++i) {
            // ALWAYS j = i+1 (AoS ile birebir)
            for (std::size_t j = i + 1; j < N; ++j) {

                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];

                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 < 1e-24) r2 = 1e-24;

                double r = std::sqrt(r2);
                double inv_r = 1.0 / r;

                double p1m = mass[i];
                double p2m = mass[j];

                const double A = _A;
                const double a1 = _a1;
                const double a2 = _a2;
                const double a_m1 = _a_m1;
                const double b = _b;
                const double C6 = _C6;
                const double C8 = _C8;
                const double C10 = _C10;

const double x0 = fast_pow(C8, 15);
        const double x1 = b*r;
        const double x2 = std::exp(x1);
        const double x3 = fast_pow(C10, 10)*fast_pow(C6, 6);
        const double x4 = fast_pow(C10, 3)*C6*fast_pow(C8, 12);
        const double x5 = fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7);
        const double x6 = fast_pow(r, 2);
        const double x7 = A*std::exp(a1*r + a2*x6 + a_m1*inv_r);
        const double x8 = fast_pow(inv_r, 2);
        const double x9 = fast_pow(b, 13);
        const double x10 = fast_pow(r, 12);
        const double x11 = b + 11*inv_r;
        const double x12 = b + inv_r;
        const double x13 = fast_pow(b, 11);
        const double x14 = fast_pow(r, 11)*x13;
        const double x15 = b + 10*inv_r;
        const double x16 = fast_pow(b, 2)*x6;
        const double x17 = b + 9*inv_r;
        const double x18 = fast_pow(b, 3)*fast_pow(r, 3);
        const double x19 = b + 8*inv_r;
        const double x20 = fast_pow(b, 4)*fast_pow(r, 4);
        const double x21 = b + 7*inv_r;
        const double x22 = fast_pow(b, 5)*fast_pow(r, 5);
        const double x23 = b + 6*inv_r;
        const double x24 = fast_pow(r, 6);
        const double x25 = fast_pow(b, 6)*x24;
        const double x26 = b + 5*inv_r;
        const double x27 = fast_pow(b, 7);
        const double x28 = fast_pow(r, 7)*x27;
        const double x29 = b + 4*inv_r;
        const double x30 = fast_pow(r, 8);
        const double x31 = fast_pow(b, 8)*x30;
        const double x32 = b + 3*inv_r;
        const double x33 = fast_pow(b, 9);
        const double x34 = fast_pow(r, 9)*x33;
        const double x35 = b + 2*inv_r;
        const double x36 = fast_pow(r, 10);
        const double x37 = fast_pow(b, 10)*x36;
        const double x38 = fast_pow(b, 15);
        const double x39 = fast_pow(r, 14);
        const double x40 = b + 13*inv_r;
        const double x41 = fast_pow(r, 13)*x9;
        const double x42 = b + 12*inv_r;
        const double x43 = fast_pow(b, 12)*x10;
        const double x44 = 720*b;
        const double x45 = 40320*b;
        const double x46 = 3628800*b;
        const double exp_tmp_0 = std::exp(-x1);

                double Fmag = -1.0/20922789888000.0*exp_tmp_0*(334764638208000*fast_pow(inv_r, 17)*x2*x3 - fast_pow(inv_r, 16)*x3*(fast_pow(b, 17)*fast_pow(r, 16) + 240*fast_pow(b, 14)*x35*x39 + 20922789888000*b + 334764638208000*inv_r + 16*fast_pow(r, 15)*x12*x38 + 20922789888000*x1*(b + 15*inv_r) + 174356582400*x11*x22 + 524160*x14*x26 + 29059430400*x15*x25 + 10461394944000*x16*(b + 14*inv_r) + 4151347200*x17*x28 + 3487131648000*x18*x40 + 518918400*x19*x31 + 871782912000*x20*x42 + 57657600*x21*x34 + 5765760*x23*x37 + 43680*x29*x43 + 3360*x32*x41) + 292919058432000*fast_pow(inv_r, 15)*x2*x5 - 240*fast_pow(inv_r, 14)*x5*(87178291200*b + 1220496076800*inv_r + 87178291200*x1*x40 + 14529715200*x11*x18 + 14*x12*x41 + 2184*x14*x32 + 3632428800*x15*x20 + 43589145600*x16*x42 + 726485760*x17*x22 + 121080960*x19*x25 + 17297280*x21*x28 + 2162160*x23*x31 + 240240*x26*x34 + 24024*x29*x37 + 182*x35*x43 + x38*x39) + 251073478656000*fast_pow(inv_r, 13)*x2*x4 - 43680*fast_pow(inv_r, 12)*x4*(479001600*b + 5748019200*inv_r + 479001600*x1*x11 + x10*x9 + 12*x12*x14 + 239500800*x15*x16 + 79833600*x17*x18 + 19958400*x19*x20 + 3991680*x21*x22 + 665280*x23*x25 + 95040*x26*x28 + 11880*x29*x31 + 1320*x32*x34 + 132*x35*x37) - 5765760*fast_pow(inv_r, 6)*x0*(C10*fast_pow(inv_r, 4)*(36288000*inv_r + r*x17*x46 + 10*x12*x34 + x13*x36 + 1814400*x16*x19 + 604800*x18*x21 + 151200*x20*x23 + 30240*x22*x26 + 5040*x25*x29 + 720*x28*x32 + 90*x31*x35 + x46) + 5040*C6*(4320*inv_r + r*x26*x44 + 6*x12*x22 + 360*x16*x29 + 120*x18*x32 + 30*x20*x35 + x24*x27 + x44) + 90*C8*x8*(322560*inv_r + r*x21*x45 + 8*x12*x28 + 20160*x16*x23 + 6720*x18*x26 + 1680*x20*x29 + 336*x22*x32 + 56*x25*x35 + x30*x33 + x45)) + 20922789888000*x0*x2*(10*C10*fast_pow(inv_r, 11) + 6*C6*fast_pow(inv_r, 7) + 8*C8*fast_pow(inv_r, 9) + a1*x7 + 2*a2*r*x7 - a_m1*x7*x8))/x0;

                double Fx = Fmag * dx * inv_r;
                double Fy = Fmag * dy * inv_r;
                double Fz = Fmag * dz * inv_r;

                fx[i] += Fx; fy[i] += Fy; fz[i] += Fz;

                if (_newton3) {
                    fx[j] -= Fx;
                    fy[j] -= Fy;
                    fz[j] -= Fz;
                }
            }
        }
    }

private:
    bool _newton3;
    double _A;
    double _a1;
    double _a2;
    double _a_m1;
    double _b;
    double _C6;
    double _C8;
    double _C10;
};
