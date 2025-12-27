def generate_soa(cfg: dict, out_dir: str):

    classname = cfg["classname"]
    params = cfg["parameters"]
    param_names = list(params.keys())

    # Force expression
    add_dispersion = all(k in param_names for k in ["C6", "C8", "C10"])

    temps_code, force_expr = calculate_force(cfg["expr_str"], param_names, add_dispersion)

    # Load templates
    header_tpl = lookup.get_template("soa_functor.h.mako")
    uses_mass = ("p1m" in cfg["expr_str"]) or ("p2m" in cfg["expr_str"])


    # Context for templates
    ctx = {
        "classname": classname,
        "parameters": params,
        "temps_code": temps_code,
        "force_expr": force_expr,
        "newton3_default": cfg["newton3_default"],
        "eps_guard": cfg["eps_guard"],
        "uses_mass": uses_mass, 
    }

    # Render
    header_code = header_tpl.render(**ctx)
    

    # Output
    hpath = os.path.join(out_dir, f"{classname}_SoA.h")
   

    with open(hpath, "w") as f:
        f.write(header_code)
    

    print("[SoA] Generated:", hpath)
