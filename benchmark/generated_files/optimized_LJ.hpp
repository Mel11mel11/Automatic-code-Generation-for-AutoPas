#pragma once
namespace lj {
    // Magnitude of Lennard-Jones force WITHOUT std::pow, keeping r-based (sqrt stays at caller).
    // Parameters are expected to be ctor-hoisted: eps24 = 24*epsilon, sigma6 = sigma^6.
    double computeForceMag_noPow(double r, double eps24, double sigma6);
} // namespace lj
