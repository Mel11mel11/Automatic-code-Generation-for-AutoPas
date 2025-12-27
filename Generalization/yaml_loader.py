import yaml
import os
import sys
# this function reads the test_yaml data and returns a list of potential definitions
def load_yaml(file_path):
    try:
        with open(file_path, "r") as f:
            docs = list(yaml.safe_load_all(f))
    except OSError as e: # file not found or cannot be opened
        raise RuntimeError(f"Cannot open YAML file: {file_path}") from e

    items = []  # empty list to store potential definitions

    for doc in docs: # iterate through each document in the YAML file
        if not doc:
            continue

        if isinstance(doc, dict) and "potentials" in doc:
            for p in doc["potentials"]: # iterate through each potential definition
                items.append(p) # add potential definition to the list

        elif isinstance(doc, dict): # single potential definition as a dictionary
            items.append(doc) # add potential definition to the list

    if not items:
        raise ValueError( 
            f"No valid potential definitions found in YAML file: {file_path}" # raise error if no definitions found
        )

    return items
