#pragma once
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include "potentials/generated_force.hpp"
#include "potentials/generated_gravity.hpp"

// This test file validates automatically generated force functions for Lennard-Jones and gravitational potentials.
//  work in reduced (dimensionless) units first,
// then hochskalieren back to SI units when reporting.
//  Small assertion utilities for test checking
// Throw an error if a condition is false
inline void assert_true(bool cond, const char* msg) {
    if (!cond) throw std::runtime_error(msg);
}
// Compare two floating point values with absolute and relative tolerances
inline void assert_near(double got, double exp,
                        double abs_tol=1e-12, double rel_tol=1e-9,
                        const char* msg="assert_near failed") {
    double diff = std::fabs(got - exp);
    double tol = std::max(abs_tol, rel_tol * std::max(1.0, std::fabs(exp)));
    if (diff > tol) throw std::runtime_error(msg);
}
// Simple check to ensure the number is finite (not NaN or infinity)
inline bool isfinite_all(double x) {
    return std::isfinite(x) && !std::isnan(x);
}
//  Reference (reduced) LJ force 
// Analytical formula in reduced units: F* = 24 * (2/r*^13 - 1/r*^7)
inline double refLJ_reduced(double rstar) {
    double inv = 1.0 / rstar;
    double inv2 = inv*inv;
    double inv6 = inv2*inv2*inv2;
    double inv12 = inv6*inv6;
    return 24.0 * (2.0*inv12 - inv6) * inv;
}
// Optional helper for physical unit conversions 
// Converts between reduced quantities and SI values.
struct LJUnits {
    double sigma;   // length scale [m]
    double epsilon; // energy scale [J]
    double mass;    // mass scale [kg]

    // Characteristic time scale τ = σ * sqrt(m / ε)
    double tau() const { return sigma * std::sqrt(mass / epsilon); }

    // Convert reduced force to SI: F = (ε/σ) * F*
    double F_from_red(double Fstar) const { return (epsilon / sigma) * Fstar; }
};
// LENNARD-JONES TESTS
inline void run_LJ_tests() {
    //  1- Force must be zero at equilibrium distance r_min ---
    {
        const double eps = 1.0, sig = 1.0;           // reduced parameters
        const double rmin = std::pow(2.0, 1.0/6.0) * sig;  // theoretical minimum
        double F = lj::computeForce(rmin, eps, sig);       // generated function
        assert_near(F, 0.0, 1e-12, 0.0, "LJ: F(r_min) must be zero");
    }
    // 2- Check the sign of force around equilibrium ---
    // For r < r_min → repulsive (positive force)
    // For r > r_min → attractive (negative force)
    {
    const double eps = 1.0, sig = 1.0;
    const double rmin = std::pow(2.0, 1.0/6.0) * sig;
    double Fl = lj::computeForce(0.95*rmin, eps, sig);
    double Fr = lj::computeForce(1.05*rmin, eps, sig);
    const double tiny = 1e-14; // 10**-14
    assert_true(std::fabs(Fl) > tiny, "LJ: left-of-min force too small?");
    assert_true(std::fabs(Fr) > tiny, "LJ: right-of-min force too small?");
    assert_true(Fl > 0.0, "LJ: repulsive region should be positive");
    assert_true(Fr < 0.0, "LJ: attractive region should be negative");
}   

    // 3- Compare generated force to analytical reference (reduced form) ---
    // Checks the shape and numerical precision over a wide range.
    {
        const double eps = 1.0, sig = 1.0;
        int N = 200;
        for (int i = 0; i < N; ++i) {
            double rstar = 0.8 + 4.2 * (double)i / (N - 1);  // from 0.8σ to 5σ
            double r = rstar * sig;                          // convert to dimensional
            double F_generated = lj::computeForce(r, eps, sig);
            double F_reference = (eps / sig) * refLJ_reduced(rstar);
            assert_near(F_generated, F_reference, 1e-12, 1e-9,
                        "LJ: analytic reduced shape comparison");
            assert_true(isfinite_all(F_generated), "LJ: value must be finite");
        }
    }
// --- 4) Scaling (homogeneity) test ---
// Keep the reduced distance r* constant.
// If σ' = a·σ and ε' = b·ε and r' = a·r  →  F' = (b/a) · F
{
    double eps = 2.3, sig = 1.7;
    // pick a reduced distance r*
    double rstar = 1.37;
    double r     = rstar * sig;
    // base force
    double F1 = lj::computeForce(r, eps, sig);
    // scale parameters
    double a = 2.0, b = 3.0;
    double sig2 = a * sig;       // σ'
    double r2   = rstar * sig2;  // r' = a·r so that r* stays the same
    // scaled force
    double F2 = lj::computeForce(r2, b * eps, sig2);
    // expectation: F2 = (b/a) * F1
    assert_near(F2, (b/a) * F1, 1e-12, 1e-10, "LJ: scaling homogeneity check");
}
    // 5- Asymptotic behavior ---
    // For large r*, the LJ force should decay ~ -24 / r*^7
    {
        const double eps = 1.0, sig = 1.0;
        for (double rstar : {10.0, 20.0, 40.0}) {
            double r = rstar * sig;
            double F_generated = lj::computeForce(r, eps, sig);
            double F_asymptotic = (eps/sig) * (-24.0 * std::pow(rstar, -7.0));
            // Allow a looser tolerance since this is asymptotic
            assert_near(F_generated, F_asymptotic, 1e-14, 5e-3,
                        "LJ: asymptotic r^-7 behavior");
        }
    }
     // 6- Numerical stability near the repulsive core ---
    // Forces must remain finite and non-NaN for small r* (0.9–1.0)
    {
        const double eps = 1.0, sig = 1.0;
        for (double rstar : {0.9, 0.92, 0.95, 1.0}) {
            double r = rstar * sig;
            double F_generated = lj::computeForce(r, eps, sig);
            assert_true(isfinite_all(F_generated), "LJ: stable near core region");
        }
    }

    // 7- Optional: reduced → SI check using Argon parameters ---
    // Demonstrates how to scale back to real-world units.
    {
        const double kB = 1.380649e-23; // Boltzmann constant [J/K]
        LJUnits ar { 3.405e-10, kB * 119.8, 6.6335209e-26 }; // Argon parameters
        double rstar = 1.3;                       // reduced distance
        double Fstar_ref = refLJ_reduced(rstar);  // reduced force
        double F_ref_SI = ar.F_from_red(Fstar_ref); // convert to Newtons

        double r = rstar * ar.sigma;              // convert distance to meters
        double F_gen_SI = lj::computeForce(r, ar.epsilon, ar.sigma);

        // Check that generated SI force matches scaled reference
        assert_near(F_gen_SI, F_ref_SI, 1e-20, 1e-9, "LJ: reduced→SI scaling test");
    }

    std::cout << "LJ tests passed 🟢\n";
}
// GRAVITATIONAL TESTS
// Reference gravitational force (SI): F = -G * m1 * m2 / r^2
inline double refGrav(double r, double G, double m1, double m2) {
    return -G * m1 * m2 / (r * r);
}
inline void run_Grav_tests() {
    const double G = 6.674e-11; // gravitational constant [m^3/(kg·s^2)]
    const double m1 = 1.0, m2 = 2.0; // example masses [kg]

    // --- 1) Compare generated vs analytical formula over a wide range ---
    {
        int N = 200;
        for (int i = 0; i < N; ++i) {
            double r = 1.0 + (1e6 - 1.0) * (double)i / (N - 1); // 1 → 10^6 m
            double F_generated = grav::computeForce(r, G, m1, m2);
            double F_reference = refGrav(r, G, m1, m2);
            assert_near(F_generated, F_reference, 1e-22, 1e-12,
                        "Grav: analytic match over range");
            assert_true(isfinite_all(F_generated), "Grav: force must be finite");
        }
    }

    // --- 2) Homogeneity (scaling) checks ---
    // Scaling distance: r' = a*r → F' = F / a^2
    // Scaling mass: m' = c*m → F' = c^2 * F
    {
        double r = 1234.5;
        double F0 = grav::computeForce(r, G, m1, m2);

        double a = 3.0; // distance scale factor
        double F_rscaled = grav::computeForce(a*r, G, m1, m2);
        assert_near(F_rscaled, F0 / (a*a), 1e-22, 1e-12,
                    "Grav: distance scaling 1/r^2 law");

        double c = 4.0; // mass scale factor
        double F_mscaled = grav::computeForce(r, G, c*m1, c*m2);
        assert_near(F_mscaled, (c*c)*F0, 1e-22, 1e-12,
                    "Grav: mass scaling m1*m2 proportionality");
    }

    // --- 3) Check that the force magnitude becomes negligible for large r ---
    {
        double r_big = 1e9; // very large distance
        double Fg = grav::computeForce(r_big, G, m1, m2);
        assert_true(std::fabs(Fg) < 1e-20, "Grav: force vanishes at large distance");
    }

    std::cout << "Gravity tests passed 🟢\n";
}

inline void run_all_tests() {
    run_LJ_tests();
    run_Grav_tests();
    std::cout << "All tests OK\n";
}
