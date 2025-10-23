#pragma once
#ifndef GENERATED_MIE_HPP
#define GENERATED_MIE_HPP

#include <stdexcept>

namespace mie {
double C(double n, double m);
double rmin(double sig, double n, double m);
double force_fast(double r, double eps, double sig, double n, double m);
double force_safe(double r, double eps, double sig, double n, double m);
} // namespace mie

#endif // GENERATED_MIE_HPP
