import os

from rdkit import Chem
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob

from hfrpkg.common import conformer_name, resolve_software, write_index
from hfrpkg.conformers import ConformerSettings, generate_conformers
from hfrpkg.core import Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic
from hfrpkg.utils import display_reaction_counts, geom_from_rdkit

reaction_map = {
    "homodesmotic": Homodesmotic,
    "isodesmic": Isodesmic,
    "isogyric": Isogyric,
    "hypohomodesmotic": Hypohomodesmotic,
}


def _smiles(m):
    return Chem.MolToSmiles(Chem.RemoveHs(m))


def _input_last(lhs, input_inchi):
    """compute_folder identifies the input by InChI, but keeping it as the last
    reactant also keeps legacy folders and eyeballing the index sane."""
    others = [(m, c) for m, c in lhs if Chem.MolToInchi(m) != input_inchi]
    mine = [(m, c) for m, c in lhs if Chem.MolToInchi(m) == input_inchi]
    return others + mine


def run_reaction(action_type, reaction_type, input_smiles, lhs=None, rhs=None, substruct=None,
                 replacement=None, outfolder=None, method=None, basis=None, extension=None, conformers=None):
    parsed = Chem.MolFromSmiles(input_smiles)
    if parsed is None:
        raise ValueError(f"Invalid SMILES: {input_smiles}")
    mol = Chem.AddHs(parsed)
    input_inchi = Chem.MolToInchi(mol)

    software, in_ext, _ = resolve_software(extension)
    level = Theory(
        method=method or "B3LYP",
        basis=basis or "6-31G",
        job_type=[OptimizationJob(), FrequencyJob()],
    )

    reaction_fn = reaction_map[reaction_type.lower()]
    rhs_mols, lhs_mols, status = reaction_fn(mol, lhs, rhs, substruct, replacement)
    if status != "Optimal":
        print(f"No balanced reaction found (solver status: {status})")
        return None
    lhs_mols = _input_last(lhs_mols, input_inchi)

    if action_type == "count":
        display_reaction_counts(mol, reaction_fn)

    elif action_type == "view":
        print("-----------Reactants-----------")
        for m, coeff in lhs_mols:
            print(f"({coeff})*{_smiles(m)}")
        print("-----------Products-----------")
        for m, coeff in rhs_mols:
            print(f"({coeff})*{_smiles(m)}")

    elif action_type == "write":
        if not outfolder:
            raise ValueError("Outfolder must be provided for write action")
        os.makedirs(outfolder, exist_ok=True)
        settings = conformers or ConformerSettings()
        entries, counts = [], {}
        for prefix, group in (("R", lhs_mols), ("P", rhs_mols)):
            for n, (m, coeff) in enumerate(group, start=1):
                name = f"{prefix}{n}_{coeff}"
                inchi = Chem.MolToInchi(m)
                mh, cids = generate_conformers(inchi, settings)
                if len(cids) == 1:  # legacy layout: one file, no suffix
                    geom_from_rdkit(mh, cids[0]).write(outfile=os.path.join(outfolder, name + in_ext), theory=level)
                else:
                    counts[name] = len(cids)
                    for k, cid in enumerate(cids, start=1):
                        geom_from_rdkit(mh, cid).write(
                            outfile=os.path.join(outfolder, conformer_name(name, k) + in_ext), theory=level)
                entries.append((name, inchi, _smiles(m)))
        write_index(outfolder, reaction_fn.__name__, software, input_smiles, entries,
                    conformers=counts, conformer_rule=settings.describe())
        print("Reaction written to", outfolder,
              f"({sum(counts.values()) or len(entries)} input files; conformers: {settings.describe()})")

    print("Complete")
    return rhs_mols, lhs_mols
