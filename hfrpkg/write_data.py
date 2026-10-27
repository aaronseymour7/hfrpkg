import os

from hfrpkg.common import KJ_PER_KCAL


def kjTokcal(value):
    if value is None or value == "":
        return None
    try:
        return value / KJ_PER_KCAL
    except TypeError:
        return None


def write_summary(reaction_data, folder_path, filename="opt_summary.txt"):
    """Write the plain-text summary and return the CSV row for it."""
    output_path = os.path.join(folder_path, filename)
    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write(f"{reaction_data['input_smiles']}\t{reaction_data['input_inchi']}\n")
        fout.write(f"Reaction Enthalpy (kJ/mol):{reaction_data['reaction_H']}\n")
        fout.write(f"DFT Enthalpy of Formation (kJ/mol):{reaction_data['dft_hf']}\n")
        for title, key in (("REACTANTS", "reactants"), ("PRODUCTS", "products")):
            fout.write(title + "\n")
            for coeff, smiles, inchi, atct, energy, zpve in reaction_data[key]:
                fout.write(f"{coeff} {smiles}\t{inchi}\t{atct}\t{energy}\t{zpve}\n")
        fout.write("\n")

    return (
        reaction_data["input_smiles"],
        reaction_data["level"],
        reaction_data["input_inchi"],
        kjTokcal(reaction_data["dft_hf"]),
        kjTokcal(reaction_data["input_atct"]),
    )


def write_single_reaction(reaction_data, folder_path):
    return write_summary(reaction_data, folder_path, "opt_summary.txt")


def write_conformer_selection(reaction_data, folder_path):
    """conformer_selection.csv: the energy of every conformer and the one that was used."""
    sel = reaction_data.get("conformers") or {}
    if not sel:
        return None
    path = os.path.join(folder_path, "conformer_selection.csv")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("species,conformer,E_plus_ZPVE_hartree,selected\n")
        for name, info in sel.items():
            for stem, e in info["energies"]:
                fh.write(f"{name},{stem},{e:.8f},{'yes' if stem == info['selected'] else ''}\n")
    return path


CSV_HEADER = ["SMILES", "LEVEL", "InChI", "ΔHf DFT (kcal/mol)", "ΔHf ATcT (kcal/mol)"]
