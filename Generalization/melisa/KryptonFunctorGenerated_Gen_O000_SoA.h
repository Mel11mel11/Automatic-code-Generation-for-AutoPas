#pragma once
#include "../Particle.h"
#include <cmath>
#include <cstddef>
#include <algorithm>  // std::max

template <class SoAView_T>
class KryptonFunctorGenerated_Gen_O000_SoA {
public:
    using SoAView = SoAView_T;

    explicit KryptonFunctorGenerated_Gen_O000_SoA(
        double A,
        double a1,
        double a2,
        double a_m1,
        double b,
        double C6,
        double C8,
        double C10,
        bool newton3 = true,
        double cutoff = 0.0
    )
      : _newton3(newton3)
      , _cutoff(cutoff)
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

        // --- SoA raw pointers ---
        double * __restrict__ x  = soa.x.data();
        double * __restrict__ y  = soa.y.data();
        double * __restrict__ z  = soa.z.data();

        double * __restrict__ fx = soa.fx.data();
        double * __restrict__ fy = soa.fy.data();
        double * __restrict__ fz = soa.fz.data();


        const std::size_t N = soa.size();
        constexpr double EPS = 1e-24;

        const double cutoff  = _cutoff;
        const double cutoff2 = cutoff * cutoff;

        const double A = _A;
        const double a1 = _a1;
        const double a2 = _a2;
        const double a_m1 = _a_m1;
        const double b = _b;
        const double C6 = _C6;
        const double C8 = _C8;
        const double C10 = _C10;

        for (std::size_t i = 0; i < N; ++i) {

            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;

            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];


            if (_newton3) {
                // =====================================================
                // Newton3 = ON  → half-shell (i < j), symmetric update
                // =====================================================
                for (std::size_t j = i + 1; j < N; ++j) {

                    const double dx = xi - x[j];
                    const double dy = yi - y[j];
                    const double dz = zi - z[j];

                    double r2 = dx*dx + dy*dy + dz*dz;
                    r2 = std::max(r2, EPS);

                    if (cutoff > 0.0 && r2 > cutoff2) continue;

                    const double r = std::sqrt(r2);
                    const double inv_r = 1.0 / r;


                    const double tt_x = b*r;
        const double exp_m_tt_x = std::exp(-tt_x);
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

                    const double Fmag = -A*(a1 + 2*a2*r - a_m1/fast_pow(r, 2))*std::exp(a1*r + a2*fast_pow(r, 2) + a_m1/r) + fast_pow(C10, 10)*fast_pow(C6, 6)*(-b*tt15*std::exp(-b*r) + b*tt16*std::exp(-b*r))/(fast_pow(C8, 15)*fast_pow(r, 16)) - 16*fast_pow(C10, 10)*fast_pow(C6, 6)*(-tt16*std::exp(-b*r) + 1)/(fast_pow(C8, 15)*fast_pow(r, 17)) + fast_pow(C10, 6)*fast_pow(C6, 3)*(-b*tt13*std::exp(-b*r) + b*tt14*std::exp(-b*r))/(fast_pow(C8, 8)*fast_pow(r, 14)) - 14*fast_pow(C10, 6)*fast_pow(C6, 3)*(-tt14*std::exp(-b*r) + 1)/(fast_pow(C8, 8)*fast_pow(r, 15)) + fast_pow(C10, 3)*C6*(-b*tt11*std::exp(-b*r) + b*tt12*std::exp(-b*r))/(fast_pow(C8, 3)*fast_pow(r, 12)) - 12*fast_pow(C10, 3)*C6*(-tt12*std::exp(-b*r) + 1)/(fast_pow(C8, 3)*fast_pow(r, 13)) + C10*(b*tt10*std::exp(-b*r) - b*tt9*std::exp(-b*r))/fast_pow(r, 10) - 10*C10*(-tt10*std::exp(-b*r) + 1)/fast_pow(r, 11) + C6*(-b*tt5*std::exp(-b*r) + b*tt6*std::exp(-b*r))/fast_pow(r, 6) - 6*C6*(-tt6*std::exp(-b*r) + 1)/fast_pow(r, 7) + C8*(-b*tt7*std::exp(-b*r) + b*tt8*std::exp(-b*r))/fast_pow(r, 8) - 8*C8*(-tt8*std::exp(-b*r) + 1)/fast_pow(r, 9);

                    const double Fx = Fmag * dx * inv_r;
                    const double Fy = Fmag * dy * inv_r;
                    const double Fz = Fmag * dz * inv_r;

                    fix += Fx;
                    fiy += Fy;
                    fiz += Fz;

                    fx[j] -= Fx;
                    fy[j] -= Fy;
                    fz[j] -= Fz;
                }

            } else {
                // =====================================================
                // Newton3 = OFF → full-shell (i != j), single-sided
                // =====================================================
                for (std::size_t j = 0; j < N; ++j) {
                    if (j == i) continue;

                    const double dx = xi - x[j];
                    const double dy = yi - y[j];
                    const double dz = zi - z[j];

                    double r2 = dx*dx + dy*dy + dz*dz;
                    r2 = std::max(r2, EPS);

                    if (cutoff > 0.0 && r2 > cutoff2) continue;

                    const double r = std::sqrt(r2);
                    const double inv_r = 1.0 / r;


                    const double tt_x = b*r;
        const double exp_m_tt_x = std::exp(-tt_x);
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

                    const double Fmag = -A*(a1 + 2*a2*r - a_m1/fast_pow(r, 2))*std::exp(a1*r + a2*fast_pow(r, 2) + a_m1/r) + fast_pow(C10, 10)*fast_pow(C6, 6)*(-b*tt15*std::exp(-b*r) + b*tt16*std::exp(-b*r))/(fast_pow(C8, 15)*fast_pow(r, 16)) - 16*fast_pow(C10, 10)*fast_pow(C6, 6)*(-tt16*std::exp(-b*r) + 1)/(fast_pow(C8, 15)*fast_pow(r, 17)) + fast_pow(C10, 6)*fast_pow(C6, 3)*(-b*tt13*std::exp(-b*r) + b*tt14*std::exp(-b*r))/(fast_pow(C8, 8)*fast_pow(r, 14)) - 14*fast_pow(C10, 6)*fast_pow(C6, 3)*(-tt14*std::exp(-b*r) + 1)/(fast_pow(C8, 8)*fast_pow(r, 15)) + fast_pow(C10, 3)*C6*(-b*tt11*std::exp(-b*r) + b*tt12*std::exp(-b*r))/(fast_pow(C8, 3)*fast_pow(r, 12)) - 12*fast_pow(C10, 3)*C6*(-tt12*std::exp(-b*r) + 1)/(fast_pow(C8, 3)*fast_pow(r, 13)) + C10*(b*tt10*std::exp(-b*r) - b*tt9*std::exp(-b*r))/fast_pow(r, 10) - 10*C10*(-tt10*std::exp(-b*r) + 1)/fast_pow(r, 11) + C6*(-b*tt5*std::exp(-b*r) + b*tt6*std::exp(-b*r))/fast_pow(r, 6) - 6*C6*(-tt6*std::exp(-b*r) + 1)/fast_pow(r, 7) + C8*(-b*tt7*std::exp(-b*r) + b*tt8*std::exp(-b*r))/fast_pow(r, 8) - 8*C8*(-tt8*std::exp(-b*r) + 1)/fast_pow(r, 9);

                    fix += Fmag * dx * inv_r;
                    fiy += Fmag * dy * inv_r;
                    fiz += Fmag * dz * inv_r;
                }
            }

            fx[i] += fix;
            fy[i] += fiy;
            fz[i] += fiz;
        }
    }

private:
    bool _newton3;
    double _cutoff;

    double _A;
    double _a1;
    double _a2;
    double _a_m1;
    double _b;
    double _C6;
    double _C8;
    double _C10;
};
