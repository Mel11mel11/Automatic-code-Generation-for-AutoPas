#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max

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


    void SoAFunctor(SoAView &soa) {

        // --- RESTRICT pointers (SIMD-friendly) ---
        double * __restrict__ x  = soa.x.data();
        double * __restrict__ y  = soa.y.data();
        double * __restrict__ z  = soa.z.data();

        double * __restrict__ fx = soa.fx.data();
        double * __restrict__ fy = soa.fy.data();
        double * __restrict__ fz = soa.fz.data();


        const std::size_t N = soa.size();
        constexpr double EPS = 1e-24;

        // --- parameter aliases (loop-invariant) ---
        const double A = _A;
        const double a1 = _a1;
        const double a2 = _a2;
        const double a_m1 = _a_m1;
        const double b = _b;
        const double C6 = _C6;
        const double C8 = _C8;
        const double C10 = _C10;

        for (std::size_t i = 0; i < N; ++i) {

            // --- local force accumulator for particle i ---
            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            // --- cache coordinates of i ---
            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];


            // Newton3 true -> j starts at i+1
            // Newton3 false -> we still start at i+1 but we must update both i and j
            for (std::size_t j = i + 1; j < N; ++j) {

                // cache coordinates of j (helps memory subsystem a bit)
                const double xj = x[j];
                const double yj = y[j];
                const double zj = z[j];

                const double dx = xi - xj;
                const double dy = yi - yj;
                const double dz = zi - zj;

                double r2 = dx * dx + dy * dy + dz * dz;
                // branchless guard (SIMD-friendly)
                r2 = std::max(r2, EPS);

                // cheaper: compute inv_r directly; avoid separate r variable
                const double inv_r = 1.0 / std::sqrt(r2);


                const double x0 = fast_pow(C8, 15);
        const double x1 = b/inv_r;
        const double x2 = std::exp(x1);
        const double x3 = fast_pow(C10, 10)*fast_pow(C6, 6);
        const double x4 = fast_pow(C10, 3)*C6*fast_pow(C8, 12);
        const double x5 = fast_pow(C10, 6)*fast_pow(C6, 3)*fast_pow(C8, 7);
        const double x6 = std::pow(inv_r, -2);
        const double x7 = A*std::exp((a1 + inv_r*(a2*x6 + a_m1*inv_r))/inv_r);
        const double x8 = fast_pow(inv_r, 2);
        const double x9 = fast_pow(b, 13);
        const double x10 = std::pow(inv_r, -12);
        const double x11 = b + 11*inv_r;
        const double x12 = b + inv_r;
        const double x13 = fast_pow(b, 11);
        const double x14 = x13/fast_pow(inv_r, 11);
        const double x15 = b + 10*inv_r;
        const double x16 = fast_pow(b, 2)*x6;
        const double x17 = b + 9*inv_r;
        const double x18 = fast_pow(b, 3)/fast_pow(inv_r, 3);
        const double x19 = b + 8*inv_r;
        const double x20 = fast_pow(b, 4)/fast_pow(inv_r, 4);
        const double x21 = b + 7*inv_r;
        const double x22 = fast_pow(b, 5)/fast_pow(inv_r, 5);
        const double x23 = b + 6*inv_r;
        const double x24 = std::pow(inv_r, -6);
        const double x25 = fast_pow(b, 6)*x24;
        const double x26 = b + 5*inv_r;
        const double x27 = fast_pow(b, 7);
        const double x28 = x27/fast_pow(inv_r, 7);
        const double x29 = b + 4*inv_r;
        const double x30 = std::pow(inv_r, -8);
        const double x31 = fast_pow(b, 8)*x30;
        const double x32 = b + 3*inv_r;
        const double x33 = fast_pow(b, 9);
        const double x34 = x33/fast_pow(inv_r, 9);
        const double x35 = b + 2*inv_r;
        const double x36 = std::pow(inv_r, -10);
        const double x37 = fast_pow(b, 10)*x36;
        const double x38 = fast_pow(b, 15);
        const double x39 = std::pow(inv_r, -14);
        const double x40 = b + 13*inv_r;
        const double x41 = x9/fast_pow(inv_r, 13);
        const double x42 = b + 12*inv_r;
        const double x43 = fast_pow(b, 12)*x10;
        const double x44 = 720*b;
        const double x45 = 40320*b;
        const double x46 = 3628800*b;
        const double x6 = fast_pow(inv_r, 7);
        const double x47 = fast_pow(inv_r, 9);
        const double x48 = 87178291200*b;
        const double x49 = 726485760*x22;
        const double x50 = 2184*x14;
        const double x51 = 121080960*x25;
        const double x52 = 17297280*x28;
        const double x53 = 14529715200*x18;
        const double x54 = 2162160*x31;
        const double x55 = 3632428800*x20;
        const double x56 = 240240*x34;
        const double x57 = 24024*x37;
        const double x58 = 182*x43;
        const double x59 = 14*x41;
        const double x60 = 87178291200*x1;
        const double x61 = 43589145600*x16;
        const double exp_tmp_0 = std::exp(-x1);

                const double Fmag = -1.0/20922789888000.0*exp_tmp_0*(240*fast_pow(inv_r, 6)*(1394852659200*fast_pow(inv_r, 12)*x2*x3 + 1220496076800*fast_pow(inv_r, 10)*x2*x5 + 1046139494400*fast_pow(inv_r, 8)*x2*x4 - 24024*x0*(C10*fast_pow(inv_r, 4)*(inv_r*(36288000*inv_r + 10*x12*x34 + x13*x36 + 1814400*x16*x19 + 604800*x18*x21 + 151200*x20*x23 + 30240*x22*x26 + 5040*x25*x29 + 720*x28*x32 + 90*x31*x35 + x46) + x17*x46) + 5040*C6*(inv_r*(4320*inv_r + 6*x12*x22 + 360*x16*x29 + 120*x18*x32 + 30*x20*x35 + x24*x27 + x44) + x26*x44) + 90*C8*x8*(inv_r*(322560*inv_r + 8*x12*x28 + 20160*x16*x23 + 6720*x18*x26 + 1680*x20*x29 + 336*x22*x32 + 56*x25*x35 + x30*x33 + x45) + x21*x45)) - 182*x4*x6*(479001600*b + 5748019200*inv_r + 479001600*x1*x11 + x10*x9 + 12*x12*x14 + 239500800*x15*x16 + 79833600*x17*x18 + 19958400*x19*x20 + 3991680*x21*x22 + 665280*x23*x25 + 95040*x26*x28 + 11880*x29*x31 + 1320*x32*x34 + 132*x35*x37) - x47*x5*(1220496076800*inv_r + x11*x53 + x12*x59 + x15*x55 + x17*x49 + x19*x51 + x21*x52 + x23*x54 + x26*x56 + x29*x57 + x32*x50 + x35*x58 + x38*x39 + x40*x60 + x42*x61 + x48)) - inv_r*x3*(fast_pow(b, 17) + 240*fast_pow(inv_r, 16)*(fast_pow(b, 14)*x35*x39 + 1394852659200*inv_r + x11*x49 + x15*x51 + x17*x52 + x19*x54 + x21*x56 + x23*x57 + x26*x50 + x29*x58 + x32*x59 + x40*x53 + x42*x55 + x48 + x60*(b + 15*inv_r) + x61*(b + 14*inv_r)) + 16*inv_r*x12*x38) + 20922789888000*x0*x2*(2*a2*x7 + inv_r*(10*C10*fast_pow(inv_r, 11) + 6*C6*x6 + 8*C8*x47 + a1*x7 - a_m1*x7*x8)))/(inv_r*x0);

                // Keep this form to match your symbolic expression usage (dx * inv_r)
                const double Fx = Fmag * dx * inv_r;
                const double Fy = Fmag * dy * inv_r;
                const double Fz = Fmag * dz * inv_r;

                // accumulate force on i locally
                fix += Fx;
                fiy += Fy;
                fiz += Fz;

                // update force on j:
                // - Newton3: subtract (pair computed once)
                // - no Newton3: we still need to add the opposite contribution to j
                fx[j] -= Fx;
                fy[j] -= Fy;
                fz[j] -= Fz;
            }

            // write back force for i ONCE
            fx[i] += fix;
            fy[i] += fiy;
            fz[i] += fiz;
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
