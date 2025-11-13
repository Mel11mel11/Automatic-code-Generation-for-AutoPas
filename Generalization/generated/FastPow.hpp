// FastPow.hpp
#pragma once
#include <cmath>

inline double fast_pow(double a, int e) noexcept {
    if (e == 0) return 1.0;
    if (e == 1) return a;

    // already checked in code
   /* if (e < 0) {
        const int p = -e;
        if (p == 1) return 1.0 / a;
        if (p == 2) return 1.0 / (a * a);
        if (p == 3) return 1.0 / (a * a * a);
        if (p == 4) { double t = a*a; return 1.0 / (t*t); }
        if (p == 5) { double t = a*a; return 1.0 / (t*t*a); }
        if (p == 6) { double t = a*a*a; return 1.0 / (t*t); }
        if (p == 7) { double t = a*a*a; return 1.0 / (t*t*a); }
        if (p == 8) { double t = a*a; double q = t*t; return 1.0 / (q*q); }
        if (p == 9) { double t = a*a; double q = t*t; return 1.0 / (q*q*a); }
        if (p == 10){ double t = a*a*a*a*a; return 1.0 / (t*t); }
        if (p == 11){ double t = a*a*a*a*a; return 1.0 / (t*t*a); }
        if (p == 12){ double t = a*a*a; return 1.0 / (t*t*t*t); }
    
        return std::pow(a, e);
    }*/

    //
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

    // büyük/kesirli üs -> standart
    return std::pow(a, e);
}
