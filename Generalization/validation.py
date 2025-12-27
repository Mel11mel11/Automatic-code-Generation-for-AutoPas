
def validate_one(data: dict) -> dict:
    name = data.get("name")
    expr = data.get("expression")
    out = data.get("output", {}) or {}
    classname = out.get("classname")
    filename = out.get("filename")
    options = data.get("options", {}) or {}
    params = data.get("parameters", {}) or {}

    if not all([name, expr, classname, filename]):
        raise ValueError("Missing required fields")

    generate_soa_flag = bool(options.get("soa", False))

    return {
        "name": name,
        "expr_str": expr,
        "classname": classname,
        "filename": filename,
        "newton3_default": bool(options.get("newton3", True)),
        "eps_guard": float(options.get("avoid_r2_zero", 1e-24)),
        "parameters": params,
        "generate_soa": generate_soa_flag,
    }
