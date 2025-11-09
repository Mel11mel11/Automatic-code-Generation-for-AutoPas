#!/usr/bin/env python3
"""
AutoFunGen: Automatic C++ Functor Generator for Molecular Potentials
--------------------------------------------------------------------
Reads a YAML file with potential definitions (Lennard-Jones, Mie, Gravity, etc.)
and generates optimized C++ functor classes for each.

If the potential is recognized (LJ, Mie, Gravity, r^2), an optimized
hand-tuned kernel is generated (no pow/sqrt overhead).

If the potential is unknown, it automatically falls back to symbolic
derivation using SymPy, directly printing a valid C++ expression
(with pow/sqrt preserved).
"""

import re
import yaml
from pathlib import Path
import sympy as sp
from sympy import cse, numbered_symbols
from sympy.printing.cxx import CXX11CodePrinter

# ==========================================================
#                 TEMPLATE FOR GENERATED FUNCTOR
# ==========================================================
# {classname}, {ctor_args}, {mask_block}, and {core} will be replaced
# with generated content for each potential.
TEMPLATE = """#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>

// ===========================================================
//  AUTO-GENERATED FUNCTOR CLASS (DO NOT EDIT MANUALLY)
// ===========================================================
template <class Particle_T>
class {classname} : public Functor<Particle_T> {{
public:
  // --- Constructor ---
  // Receives the potential parameters (sigma, epsilon, etc.) + newton3 flag.
  {classname}({ctor_args}) : _newton3(newton3){extra_inits} {{}}

  // --- Newton3 control functions ---
  bool allowsNewton3() const {{ return true; }}
  bool usesNewton3()  const {{ return _newton3; }}

  // --- Main computation function (AoS version) ---
  inline void AoSFunctor(Particle_T& a, Particle_T& b) override {{
    // Load particle positions
    const auto& ra = a.getR();
    const auto& rb = b.getR();

    // Compute pairwise distance components
    const double dx = ra[0]-rb[0];
    const double dy = ra[1]-rb[1];
    const double dz = ra[2]-rb[2];

    // Square distance (no sqrt)
    const double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 == 0.0) return;  // Avoid division by zero

    // Apply cutoff mask
    {mask_block}

    // --- Core force computation (auto-generated) ---
    {core}

    // Accumulate force on both particles
    a.addF(Fx, Fy, Fz);
    if (_newton3) b.subF(Fx, Fy, Fz);
  }}

private:
  bool _newton3;  // Newton3 flag (true if symmetric forces allowed)
{members}
}};
"""

# ==========================================================
#          OPTIMIZED FORCE CORES (MANUALLY TUNED)
# ==========================================================

def core_lj():
    """Optimized Lennard-Jones force core (no sqrt, single division)."""
    return """const double inv_r2 = 1.0 / r2;
    const double sr2 = _sigma2 * inv_r2;
    const double sr6 = sr2 * sr2 * sr2;
    const double fac = _eps24 * (2.0*sr6*sr6 - sr6) * inv_r2 * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;"""

def core_gravity():
    """Gravitational potential F = -G*m1*m2/r^3 (needs one sqrt)."""
    return """const double inv_r2 = 1.0 / r2;
    const double inv_r  = 1.0 / std::sqrt(r2);
    const double inv_r3 = inv_r2 * inv_r;
    const double C = _G * _m1 * _m2;
    const double fac = -C * inv_r3 * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;"""

def core_mie_even():
    """Mie potential (n,m both even) — no sqrt needed."""
    return """const double inv_r2 = 1.0 / r2;
    double sig_pow_n = 1.0, sig_pow_m = 1.0;
    { // Compute sigma^(n) and sigma^(m) via (sigma^2)^k
      double s2 = _sigma2;
      int hn = _n/2, hm = _m/2;
      for(int i=0;i<hn;i++) sig_pow_n *= s2;
      for(int i=0;i<hm;i++) sig_pow_m *= s2;
    }
    double inv_r_pow_n2 = 1.0, inv_r_pow_m2 = 1.0;
    {{
      int pn = (_n+2)/2, pm = (_m+2)/2;
      for(int i=0;i<pn;i++) inv_r_pow_n2 *= inv_r2;
      for(int i=0;i<pm;i++) inv_r_pow_m2 *= inv_r2;
    }}
    const double term_n = _n * sig_pow_n * inv_r_pow_n2;
    const double term_m = _m * sig_pow_m * inv_r_pow_m2;
    const double fac = _C * _epsilon * (term_n - term_m) * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;"""

def core_mie_general():
    """Mie potential (n or m odd) — needs one sqrt for inv_r."""
    return """const double inv_r2 = 1.0 / r2;
    const double inv_r  = 1.0 / std::sqrt(r2);
    double sigma_n = 1.0, sigma_m = 1.0;
    {{
      double s = std::sqrt(_sigma2); // sigma
      for(int i=0;i<_n;i++) sigma_n *= s;
      for(int i=0;i<_m;i++) sigma_m *= s;
    }}
    auto inv_r_pow = [&](int k){
      int full2 = k/2, rest1 = k%2;
      double t=1.0;
      for(int i=0;i<full2;i++) t*=inv_r2;
      if(rest1) t*=inv_r;
      return t;
    };
    const double term_n = _n * sigma_n * inv_r_pow(_n+2);
    const double term_m = _m * sigma_m * inv_r_pow(_m+2);
    const double fac = _C * _epsilon * (term_n - term_m) * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;"""

def core_test_r2():
    """Simple V = r^2 potential — just linear force."""
    return """const double fac = -2.0 * mask;
    const double Fx = dx * fac;
    const double Fy = dy * fac;
    const double Fz = dz * fac;"""

# ==========================================================
#          FALLBACK: GENERIC SYMBOLIC DERIVATION
# ==========================================================

def fallback_direct_core(cfg):
    """
    Generic fallback for unknown potentials.
    - Uses SymPy to compute dV/dr symbolically.
    - Keeps pow/sqrt/exp/etc. intact (no optimizations).
    - Produces valid C++ code directly.
    """
    # Define coordinate symbols
    dx, dy, dz = sp.symbols('dx dy dz', real=True)
    r2 = dx*dx + dy*dy + dz*dz
    r  = sp.sqrt(r2)

    # Load YAML parameters as symbolic variables
    params = cfg.get("parameters", {})
    ps = {k: sp.symbols(k, real=True) for k in params.keys()}

    # Allow standard math functions in expressions
    local_dict = {**ps, "r": r, "sqrt": sp.sqrt, "exp": sp.exp,
                  "log": sp.log, "sin": sp.sin, "cos": sp.cos, "tanh": sp.tanh}

    # Parse user expression (example: "A*exp(-k*r)/r")
    V = sp.sympify(cfg["expression"], locals=local_dict)

    # Compute dV/dr and convert to force factor
    dV_dr = sp.diff(V, r)
    mask  = sp.Symbol("mask", real=True)
    fac   = -dV_dr * (1/r) * mask

    # Compute Fx, Fy, Fz components
    Fx = dx * fac
    Fy = dy * fac
    Fz = dz * fac

    # Common subexpression elimination (CSE)
    repls, (Fx_c, Fy_c, Fz_c) = cse([Fx, Fy, Fz], symbols=numbered_symbols('t'))

    # Use default C++ printer (pow/sqrt remain visible)
    printer = CXX11CodePrinter({'standard': 'c++11'})

    # Emit declaration + expression lines
    decls = [f"double {str(s)} = {printer.doprint(e)};" for s, e in repls]
    lines = decls + [
        f"double Fx = {printer.doprint(Fx_c)};",
        f"double Fy = {printer.doprint(Fy_c)};",
        f"double Fz = {printer.doprint(Fz_c)};",
    ]
    return "\n    ".join(lines)

# ==========================================================
#                HELPER FUNCTIONS
# ==========================================================
def _sanitize(name: str) -> str:
    """Sanitize a parameter name so it becomes a valid C++ identifier."""
    s = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    if s and s[0].isdigit():
        s = "_" + s
    return s

def load_yaml_all(path):
    """Load all YAML documents separated by '---'."""
    text = Path(path).read_text(encoding="utf-8")
    return [d for d in yaml.safe_load_all(text) if d]

# ==========================================================
#       MAIN DISPATCHER: chooses which core to use
# ==========================================================
def detect_and_emit(cfg):
    """Main logic: decides which core generator to use and writes C++ file."""
    expr = cfg["expression"].replace(" ", "")
    params = cfg.get("parameters",{}) or {}
    classname = cfg["output"]["classname"]
    filename  = cfg["output"]["filename"]

    # --- Cutoff mask (if rcut present) ---
    if "rcut" in params:
        mask_block = "const double _rcut2 = rcut*rcut;\n    const double mask = (r2 < _rcut2) ? 1.0 : 0.0;"
    else:
        mask_block = "const double mask = 1.0;"

    ctor_args, members, extra_inits, core = [], [], "", ""

    # --- Pattern detection for known potentials ---
    if "(sigma/r)**12" in expr and "(sigma/r)**6" in expr:
        # Lennard-Jones potential
        ctor_args += ["double sigma","double epsilon"]
        members   += ["double _sigma2 = 0.0;","double _eps24 = 0.0;"]
        extra_inits += ", _sigma2(sigma*sigma), _eps24(24.0*epsilon)"
        if "rcut" in params: ctor_args += ["double rcut"]
        core = core_lj()

    elif "/r" in expr and "sigma" not in expr and "epsilon" not in expr and "(sigma/r)" not in expr:
        # Gravity-like potential: -G*m1*m2 / r
        ctor_args += ["double G","double m1","double m2"]
        members   += ["double _G = 0.0;","double _m1 = 0.0;","double _m2 = 0.0;"]
        extra_inits += ", _G(G), _m1(m1), _m2(m2)"
        if "rcut" in params: ctor_args += ["double rcut"]
        core = core_gravity()

    elif "(sigma/r)**n" in expr and "(sigma/r)**m" in expr:
        # Mie potential
        ctor_args += ["double sigma","double epsilon","int n","int m","double C"]
        members   += ["double _sigma2 = 0.0;","int _n = 0;","int _m = 0;",
                      "double _C = 0.0;","double _epsilon = 0.0;"]
        extra_inits += ", _sigma2(sigma*sigma), _n(n), _m(m), _C(C), _epsilon(epsilon)"
        if "rcut" in params: ctor_args += ["double rcut"]
        n = int(params.get("n",12)); m = int(params.get("m",6))
        core = core_mie_even() if (n%2==0 and m%2==0) else core_mie_general()

    elif expr == "r**2":
        # Simple test potential V=r²
        core = core_test_r2()

    else:
        # --- Unknown potential: use symbolic fallback ---
        core = fallback_direct_core(cfg)

        # Automatically add YAML parameters as class members + ctor args
        ctor_params = []
        members_params = []
        extra_inits_params = []
        for pname in params.keys():
            if pname == "rcut":
                continue
            safe = _sanitize(pname)
            ctor_params.append(f"double {safe}")
            members_params.append(f"double _{safe} = 0.0;")
            extra_inits_params.append(f"_{safe}({safe})")
            # Replace parameter names in code with member variables
            core = re.sub(rf'\\b{pname}\\b', f'_{safe}', core)

        ctor_args += ctor_params
        members   += members_params
        if extra_inits and extra_inits.strip():
            if extra_inits_params:
                extra_inits += ", " + ", ".join(extra_inits_params)
        else:
            extra_inits  = (", " + ", ".join(extra_inits_params)) if extra_inits_params else ""

        if "rcut" in params:
            ctor_args += ["double rcut"]

    # --- Finalize constructor arguments ---
    ctor_args.append("bool newton3=true")
    members_str = ("\n  ".join(members)) if members else ""

    # --- Render C++ class ---
    code = TEMPLATE.format(
        classname=classname,
        ctor_args=", ".join(ctor_args),
        extra_inits=extra_inits,
        members=members_str,
        mask_block=mask_block,
        core=core
    )

    # Write output file
    Path(filename).write_text(code, encoding="utf-8")
    print(f"[ok] {classname} -> {filename}")

# ==========================================================
#                       ENTRY POINT
# ==========================================================
def main(yaml_path):
    """Reads the YAML file and generates C++ functors for each potential."""
    for cfg in load_yaml_all(yaml_path):
        detect_and_emit(cfg)

if __name__ == "__main__":
    import sys
    main(sys.argv[1])
