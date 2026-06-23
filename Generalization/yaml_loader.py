import yaml


def load_yaml(file_path: str) -> list[dict]:
    """
    Read a YAML file that may contain one or more potential specs (separated by ---).
    Returns a flat list of raw potential definition dicts.
    """
    try:
        with open(file_path, "r") as f:
            docs = list(yaml.safe_load_all(f))
    except OSError as e:
        raise RuntimeError(f"Cannot open YAML file: {file_path}") from e

    items = []

    for doc in docs:
        if not doc:
            continue

        # A file may group multiple potentials under a "potentials:" key
        if isinstance(doc, dict) and "potentials" in doc:
            for p in doc["potentials"]:
                items.append(p)
        elif isinstance(doc, dict):
            items.append(doc)

    if not items:
        raise ValueError(f"No valid potential definitions found in: {file_path}")

    return items
