"""Turn a folder of finished jobs into a reaction enthalpy and a dHf for the input.

The arithmetic lives in `build_reaction`, which takes plain numbers and so can be
unit-tested without any QM output. `compute_folder` / `compute_spec_folder` only
gather energies from log files and hand them over.
"""

import os

from hfrpkg.common import (
    HARTREE_TO_KJ,
    get_Hf,
    inchi_from_smiles,
    parse_name,
    read_index,
)


class ComputeError(Exception):
    pass


def build_reaction(index, energies, folder=".", selection=None):
    """
    index    : common.Index for the folder
    energies : {name: (E_plus_zpve_hartree, zpve_hartree)}
    Returns the summary dict. Raises ComputeError if any species has no energy:
    summing a partial set would silently give a wrong reaction enthalpy.
    """
    rows = []
    for name, inchi, smiles in index.entries:
        parsed = parse_name(name)
        if parsed is None:
            continue
        kind, pos, coeff = parsed
        rows.append((kind, pos, coeff, name, inchi, smiles))
    reactants = sorted((r for r in rows if r[0] == "R"), key=lambda r: r[1])
    products = sorted((r for r in rows if r[0] == "P"), key=lambda r: r[1])
    if not reactants:
        raise ComputeError("No reactants found in index.txt")

    missing = [r[3] for r in reactants + products if energies.get(r[3], (None, None))[0] is None]
    if missing:
        raise ComputeError("No energy for: " + ", ".join(missing))

    # The input is identified by InChI, not by file order.
    input_inchi = inchi_from_smiles(index.input_smiles) if index.input_smiles else None
    input_row = next((r for r in reactants if input_inchi and r[4] == input_inchi), None)
    if input_row is None:
        input_row = reactants[-1]  # legacy folders: input was written last
    others = [r for r in reactants if r is not input_row]

    def total(group):
        return sum(r[2] * energies[r[3]][0] for r in group)

    reaction_H = HARTREE_TO_KJ * (total(products) - total(reactants))

    def record(row):
        kind, pos, coeff, name, inchi, smiles = row
        e, z = energies[name]
        return (coeff, smiles, inchi, get_Hf(inchi), e, z)

    reactants_data = [record(r) for r in others] + [record(input_row)]
    products_data = [record(r) for r in products]

    unknown = [d for d in reactants_data[:-1] + products_data if d[3] is None]
    for d in unknown:
        print(f"[ATcT MISSING] {folder}  {d[1]} -> InChI: {d[2]}")

    input_coeff = input_row[2]
    input_hf = None
    if not unknown:
        hf_p = sum(d[0] * d[3] for d in products_data)
        hf_r = sum(d[0] * d[3] for d in reactants_data[:-1])
        input_hf = round((hf_p - hf_r - reaction_H) / input_coeff, 4)

    return {
        "input_smiles": index.input_smiles or input_row[5],
        "input_inchi": input_row[4],
        "dft_hf": input_hf,
        "level": index.level,
        "reactants": reactants_data,
        "products": products_data,
        "input_atct": get_Hf(input_row[4]),
        "reaction_H": reaction_H,
        "conformers": selection or {},
    }


def _reader(logfile):
    from AaronTools.fileIO import FileReader

    return FileReader(logfile, just_geom=False)


def read_opt_energy(logfile):
    """(E + ZPVE, ZPVE) in hartree from a finished frequency job, or (None, None)."""
    try:
        r = _reader(logfile)
        keys = r.keys()
        zpve = r["ZPVE"] if "ZPVE" in keys else None
        if "E_ZPVE" in keys:
            return r["E_ZPVE"], zpve
        if "energy" in keys and zpve is not None:
            return r["energy"] + zpve, zpve
    except Exception as exc:
        print(f"[WARN] could not read {logfile}: {exc}")
    return None, None


def read_sp_energy(logfile, zpve):
    """(E_sp + ZPVE_from_opt, ZPVE) from a single-point job."""
    if zpve is None:
        return None, None
    try:
        r = _reader(logfile)
        if "energy" in r.keys():
            return r["energy"] + zpve, zpve
    except Exception as exc:
        print(f"[WARN] could not read {logfile}: {exc}")
    return None, None


def optimized_energies(folder_path, index):
    """Lowest-energy conformer of every species: ({name: (E + ZPVE, ZPVE)}, {name: selection info})."""
    from hfrpkg.conformers import select_lowest

    _, out_ext = index.extensions
    energies, selection = {}, {}
    for name, _, _ in index.entries:
        n = index.conformers.get(name, 1)
        try:
            _, energies[name], selection[name] = select_lowest(folder_path, name, n, out_ext)
        except ValueError as exc:
            if n > 1:
                raise ComputeError(str(exc))
            energies[name] = (None, None)  # reported with the other missing energies by build_reaction
    return energies, {k: v for k, v in selection.items() if v["n"] > 1}


def compute_folder(folder_path):
    """Optimisation-level result for one *.mhfr folder; None (with a message) on failure."""
    try:
        index = read_index(folder_path)
        energies, selection = optimized_energies(folder_path, index)
        return build_reaction(index, energies, folder_path, selection)
    except Exception as exc:
        print(f"[ERROR] Folder {folder_path}: {exc}")
        return None


def compute_spec_folder(folder_path):
    """Single-point result: SP energy + ZPVE from the optimisation logs."""
    try:
        index = read_index(folder_path)
        _, opt_out = index.extensions
        spec_dir = os.path.join(folder_path, "spec")
        spec_index_path = os.path.join(spec_dir, "index.txt")
        sp_software = read_index(spec_index_path).software if os.path.exists(spec_index_path) else index.software
        from hfrpkg.common import resolve_software

        sp_out = resolve_software(sp_software)[2]
        from hfrpkg.conformers import select_lowest

        energies = {}
        for name, _, _ in index.entries:
            sp_log = os.path.join(spec_dir, name + sp_out)
            try:  # ZPVE comes from the lowest-energy conformer, the geometry the single point was run on
                _, (_, zpve), _ = select_lowest(folder_path, name, index.conformers.get(name, 1), opt_out)
            except ValueError:
                energies[name] = (None, None)
                continue
            energies[name] = read_sp_energy(sp_log, zpve) if os.path.exists(sp_log) else (None, None)
        return build_reaction(index, energies, spec_dir)
    except Exception as exc:
        print(f"[ERROR] Folder {folder_path}/spec: {exc}")
        return None
