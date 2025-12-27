// FastPow.hpp
#pragma once
#include <cmath>

inline double fast_pow(double a, int e) noexcept {
    if (e == 0) return 1.0;
    if (e == 1) return a;
    if (e == 2)  return a * a;
    if (e == 3)  return a * a * a;
    if (e == 4)  { double t = a*a; return t*t; }
    if (e == 5)  { double t = a*a; return t*t*a; }
    if (e == 6)  { double t = a*a*a; return t*t; }
    if (e == 7)  { double t = a*a*a; return t*t*a; }
    if (e == 8)  { double t = a*a; double q = t*t; return q*q; }
    if (e == 9)  { double t = a*a; double q = t*t; return q*q*a; }
    if (e == 10) { double t = a*a*a*a*a; return t*t; }
    if (e == 11) { double t = a*a*a*a*a; return t*t*a; }
    if (e == 12) { double t = a*a*a; return t*t*t*t; }
    if (e == 13) { double t = a*a*a; return t*t*t*t*a; }
    if (e == 14) { double t = a*a*a*a; return t*t*t*t; }
    if (e == 15) { double t = a*a*a*a; return t*t*t*t*a; }
    if (e == 16) { double t = a*a; double q = t*t; return q*q*q*q; }
    if (e == 17) { double t = a*a; double q = t*t; return q*q*q*q*a; }
    if (e == 18) { double t = a*a*a; double q = t*t*t; return q*q; }
    if (e == 19) { double t = a*a*a; double q = t*t*t; return q*q*a; }
    if (e == 20) { double t = a*a*a*a*a; double q = t*t; return q*q; }
    
//
    return std::pow(a, e);
}
