
import sympy as sp
from pathlib import Path

# === 1) Semboller ve potansiyel tanımı ===
r, epsilon, sigma = sp.symbols('r epsilon sigma', positive=True)
V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

# === 2) Kuvvet magnitüdü (|F|) ===
# Not: |F|(r) = -dV/dr. (V(r) merkezi, yöneyi AoS'ta r_vec / r ile verirsin.)
Fmag = -sp.diff(V, r)
Fmag = sp.simplify(Fmag)

# === 3) Ctor-hoist edilmiş sabitler: eps24 = 24*epsilon, sigma6 = sigma**6 ===
eps24, sigma6 = sp.symbols('eps24 sigma6', positive=True)
Fmag_hoisted = Fmag.subs({24*epsilon: eps24, sigma**6: sigma6})

# İsteğe bağlı ek sadeleştirme: (sigma/r)^12 = (sigma^6/r^6)^2
# (sigma/r)^6  = sigma6 / r^6
Fmag_hoisted = sp.simplify(Fmag_hoisted)

# === 4) Çıktı C++ fonksiyonunu ELDEN yaz: pow yok, sadece çarpma ===
# |F| = eps24 * ( 2*(sigma6/r^6)^2 - (sigma6/r^6) ) * (1/r)
#      = eps24 * ( 2*s6r6*s6r6 - s6r6 ) * inv_r
#  s6r6 = sigma6 * inv_r6
#
# C++ içinde:
#   double inv_r  = 1.0 / r;
#   double inv_r2 = inv_r * inv_r;
#   double inv_r6 = inv_r2 * inv_r2 * inv_r2;
#   double s6_over_r6 = sigma6 * inv_r6;
#   return eps24 * (2.0 * s6_over_r6 * s6_over_r6 - s6_over_r6) * inv_r;

header_code = r"""#pragma once
namespace lj {
    // Magnitude of Lennard-Jones force WITHOUT std::pow, keeping r-based (sqrt stays at caller).
    // Parameters are expected to be ctor-hoisted: eps24 = 24*epsilon, sigma6 = sigma^6.
    double computeForceMag_noPow(double r, double eps24, double sigma6);
} // namespace lj
"""

cpp_code = r"""#include <cmath>
#include "generated_force.hpp"

namespace lj {
inline double computeForceMag_noPow_impl(double r, double eps24, double sigma6) {
    // Guard against r=0 from caller side if needed; this routine assumes r>0.
    double inv_r  = 1.0 / r;
    double inv_r2 = inv_r * inv_r;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;

    double s6_over_r6 = sigma6 * inv_r6;                // (sigma/r)^6
    return eps24 * (2.0 * s6_over_r6 * s6_over_r6 - s6_over_r6) * inv_r;
}

double computeForceMag_noPow(double r, double eps24, double sigma6) {
    // Optional tiny clamp for stability (aynı minR2 mantığını r tabanında uygulamak istersen):
    // constexpr double minR = 1e-4;  // örnek (r^2=1e-8 ile uyumlu)
    // if (r < minR) r = minR;

    return computeForceMag_noPow_impl(r, eps24, sigma6);
}
} // namespace lj
"""

# === 5) Dosyaları yaz ===
outdir = Path("../generated_files")
outdir.mkdir(parents=True, exist_ok=True)

(outdir / "optimized_LJ.hpp").write_text(header_code, encoding="utf-8")
(outdir / "optimized_LJ.cpp").write_text(cpp_code, encoding="utf-8")

print("🟢 C++ code generated: ../generated_files/optimized_LJ.hpp/.cpp")

# === 6) AoS tarafında kullanım örneği (bilgi amaçlı; dosyaya yazmıyoruz) ===
usage_hint = r"""
// C++ tarafında ctor'da:
//   _eps24  = 24.0 * epsilon;
//   _sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
//
// AoSFunctor içinde:
//   const double mag   = lj::computeForceMag_noPow(r, _eps24, _sigma6);
//   const double inv_r = 1.0 / r;
//   std::array<double,3> F{ mag * dx * inv_r, mag * dy * inv_r, mag * dz * inv_r };
"""
print(usage_hint)
