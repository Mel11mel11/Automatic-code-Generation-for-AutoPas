import textwrap


def emit_header(
    functorname: str,
    temp_calc: str,
    force_expr: str,
    newton3_default: bool,
    eps_guard: float,
    potential_param: list[str],
    add_dispersion: bool,
    use_fast_pow: bool,
    use_mass: bool,
) -> str:
    """
    Render a complete C++ AoS functor header as a string.

    Parameters
    ----------
    functorname     : class name for the generated functor
    temp_calc       : C++ lines for intermediate variables (CSE + TT prelude)
    force_expr      : final force magnitude expression in C++
    newton3_default : default value for the newton3 flag in the constructor
    eps_guard       : minimum r^2 to avoid division by zero
    potential_param : list of parameter names (sigma, epsilon, ...)
    add_dispersion  : True for Krypton-style potentials with C6/C8/C10 terms
    use_fast_pow    : whether to include the fast_pow macro and header
    use_mass        : whether particle masses (p1m, p2m) are needed (Gravity)
    """
    potential_param = potential_param or []

    # Constructor signature: all potential parameters + newton3 flag + cutoff
    ctor_params = ", ".join(f"double {p}" for p in potential_param)
    if ctor_params:
        ctor_params += ", "
    ctor_params += f"bool newton3 = {'true' if newton3_default else 'false'}"
    ctor_params += ", double cutoff = 0.0"

    # Member initializer list
    init_elems = [f"_{p}({p})" for p in potential_param]
    if add_dispersion:
        init_elems += ["_C12(0.0)", "_C14(0.0)", "_C16(0.0)"]
    if functorname.lower().startswith("mie"):
        init_elems.append("_C(0.0)")
    init_elems += ["_newton3(newton3)", "_cutoff(cutoff)"]
    init_list = ", ".join(init_elems)

    # Local aliases so the force expression can use plain parameter names
    alias_lines = [f"const double {p} = _{p};" for p in potential_param]
    if add_dispersion:
        alias_lines += [
            "const double C12 = _C12;",
            "const double C14 = _C14;",
            "const double C16 = _C16;",
        ]
    if functorname.lower().startswith("mie"):
        alias_lines.append("const double C = _C;")
    local_aliases = "\n        ".join(alias_lines)

    # Private member declarations
    member_lines = [f"double _{p};" for p in potential_param]
    if add_dispersion:
        member_lines += ["double _C12;", "double _C14;", "double _C16;"]
    if functorname.lower().startswith("mie"):
        member_lines.append("double _C;")
    member_lines += ["bool _newton3;", "double _cutoff;"]
    member_decls = "\n    ".join(member_lines)

    # Mie potential requires computing the prefactor C at construction time
    mie_runtime = ""
    if functorname.lower().startswith("mie"):
        mie_runtime = """
        _C = (_n / (_n - _m)) * std::pow(_n / _m, _m / (_n - _m));
        """

    # Krypton dispersion: C12, C14, C16 are derived from the base C6/C8/C10 parameters
    disp_runtime = ""
    if add_dispersion:
        disp_runtime = """
        _C12 = _C6 * std::pow(_C10 / _C8, 3.0);
        _C14 = _C8 * std::pow(_C12 / _C10, 3.0);
        _C16 = _C10 * std::pow(_C14 / _C12, 3.0);
        """

    fastpow_define = "#define USE_FAST_POW\n" if use_fast_pow else ""
    mass_define = "#define USE_MASS\n" if use_mass else ""

    body = f"""
#pragma once
{fastpow_define}{mass_define}
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>

#ifdef USE_FAST_POW
#include "FastPow.hpp"
#endif

template <class Particle_T>
class {functorname} : public Functor<Particle_T> {{
public:
    explicit {functorname}({ctor_params})
        : {init_list}
    {{
        {mie_runtime}
        {disp_runtime}
    }}

    void AoSFunctor(Particle_T& p1, Particle_T& p2) override {{
        const auto& ra = p1.getR();
        const auto& rb = p2.getR();
        double dx = ra[0] - rb[0];
        double dy = ra[1] - rb[1];
        double dz = ra[2] - rb[2];

        constexpr double EPS = {eps_guard};
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < EPS) r2 = EPS;

        const double cutoff = _cutoff;
        const double cutoff2 = cutoff * cutoff;
        if (cutoff > 0.0 && r2 > cutoff2) return;

        const double r = std::sqrt(r2);
        const double inv_r = 1.0 / r;

        {local_aliases}

        #ifdef USE_MASS
        const double p1m = p1.getMass();
        const double p2m = p2.getMass();
        #endif

{("        " + temp_calc) if temp_calc else ""}

        const double Fmag = {force_expr};

        const double fx = Fmag * dx * inv_r;
        const double fy = Fmag * dy * inv_r;
        const double fz = Fmag * dz * inv_r;

        std::array<double,3> F{{fx, fy, fz}};
        p1.addF(F);
        if (_newton3) p2.subF(F);
    }}

    bool allowsNewton3() const {{ return true; }}
    bool usesNewton3() const {{ return _newton3; }}

private:
    {member_decls}
}};
"""
    return textwrap.dedent(body)
