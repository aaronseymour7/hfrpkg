import os


def _num(token):
    return None if token in ("None", "") else float(token)


def read_optsum(folder_path, filename="opt_summary.txt"):
    input_path = os.path.join(folder_path, filename)
    data = {
        "input_smiles": None,
        "input_inchi": None,
        "reaction_H": None,
        "dft_hf": None,
        "reactants": [],
        "products": [],
    }

    with open(input_path, "r", encoding="utf-8") as fin:
        lines = [line.rstrip("\n") for line in fin if line.strip()]

    first = lines[0].split("\t")
    data["input_smiles"] = first[0]
    data["input_inchi"] = first[1] if len(first) > 1 else None
    data["reaction_H"] = _num(lines[1].split(":", 1)[1])
    data["dft_hf"] = _num(lines[2].split(":", 1)[1])

    section = None
    for line in lines[3:]:
        if line == "REACTANTS":
            section = "reactants"
        elif line == "PRODUCTS":
            section = "products"
        elif section:
            parts = line.split("\t")
            coeff, smiles = parts[0].split(" ", 1)
            data[section].append(
                (int(coeff), smiles, parts[1], _num(parts[2]), _num(parts[3]), _num(parts[4]))
            )
    return data
