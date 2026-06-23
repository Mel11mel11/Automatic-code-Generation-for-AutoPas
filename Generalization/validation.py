def validate_one(data: dict) -> dict:
    """
    Validate and normalize a single potential spec from the YAML file.
    Raises ValueError if required fields are missing or values are out of range.
    """
    name = data.get("name")
    expr = data.get("expression")
    out = data.get("output", {}) or {}
    classname = out.get("classname")
    filename = out.get("filename")
    options = data.get("options", {}) or {}
    params = data.get("parameters", {}) or {}

    cutoff = options.get("cutoff", None)
    if cutoff is not None:
        cutoff = float(cutoff)
        if cutoff <= 0.0:
            raise ValueError("cutoff must be positive if specified")

    if not all([name, expr, classname, filename]):
        raise ValueError(
            f"Potential '{name}' is missing one of: expression, classname, filename"
        )

    return {
        "name": name,
        "expr_str": expr,
        "classname": classname,
        "filename": filename,
        "newton3_default": bool(options.get("newton3", True)),
        "eps_guard": float(options.get("avoid_r2_zero", 1e-24)),
        "parameters": params,
        "generate_soa": bool(options.get("soa", False)),
        "cutoff": cutoff,
    }
