import textwrap

def emit_header(functorname, temp_calc,force_expr,newton3_default,eps_guard,potential_param,
    add_dispersion,
):
    if potential_param:
        ctor_param_sig = ", ".join(f"double {p}" for p in potential_param) + ", "
    else:
        ctor_param_sig = ""

    ctor_param_sig += f"bool newton3 = {'true' if newton3_default else 'false'}"
    ctor_param_sig += ", double cutoff = 0.0"
    # init list
    potential_param = potential_param or []
    init_list_elems = [f"_{p}({p})" for p in potential_param]
    if add_dispersion:
        init_list_elems += ["_C12(0.0)", "_C14(0.0)", "_C16(0.0)"]
    if functorname.lower().startswith("mie"):
        init_list_elems.append("_C(0.0)")
    init_list_elems.append("_newton3(newton3)")
    init_list_elems.append("_cutoff(cutoff)")

    init_list = ", ".join(init_list_elems)

    # aliases
    alias_list = [f"const double {p} = _{p};" for p in potential_param]
    if add_dispersion:
        alias_list += [
            "const double C12 = _C12;",
            "const double C14 = _C14;",
            "const double C16 = _C16;",
        ]
    if functorname.lower().startswith("mie"):
        alias_list.append("const double C = _C;")
    local_aliases = "\n        ".join(alias_list)

    # member declarations
    member_decls = [f"double _{p};" for p in potential_param]
    if add_dispersion:
        member_decls += ["double _C12;", "double _C14;", "double _C16;"]
    if functorname.lower().startswith("mie"):
        member_decls.append("double _C;")
    member_decls.append("bool _newton3;")
    member_decls.append("double _cutoff;")

    member_decls_str = "\n    ".join(member_decls)

    # Mie runtime C compute
    if functorname.lower().startswith("mie"):
        mie_runtime = """
        _C = (_n / (_n - _m)) * std::pow(_n / _m, _m / (_n - _m));
        """
    else:
        mie_runtime = ""

    # Dispersion C12-C16 runtime compute
    if add_dispersion:
        disp_runtime = """
        _C12 = _C6 * std::pow(_C10 / _C8, 3.0);
        _C14 = _C8 * std::pow(_C12 / _C10, 3.0);
        _C16 = _C10 * std::pow(_C14 / _C12, 3.0);
        """
    else:
        disp_runtime = ""

    body = f"""
#pragma once
#include "../Functors/Functor.h"
#include "../Particle.h"
#include <cmath>
#include <array>
#include "FastPow.hpp"

template <class Particle_T>
class {functorname} : public Functor<Particle_T> {{
public:
    explicit {functorname}({ctor_param_sig})
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

        // Parameter aliases
        {local_aliases}

        const double p1m = p1.getMass();
        const double p2m = p2.getMass();

{("        " +  temp_calc) if  temp_calc else ""}

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
    {member_decls_str}
}};
"""
    return textwrap.dedent(body)
