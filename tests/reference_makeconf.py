"""Reference implementation: the conformer generation of makeConf.py used for the published benchmark.

Copied verbatim from the author's script, except that the final file-writing loop is replaced by a
`return` so the result can be compared with hfrpkg.conformers.generate_conformers.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign, rdMolDescriptors


def prune_conformers_by_rmsd_keep_order(mol, candidate_conf_ids, rmsd_cutoff=0.15):

    kept = []
    for cid in candidate_conf_ids:
        keep = True
        for kcid in kept:
            rmsd = rdMolAlign.GetBestRMS(mol, mol, prbId=cid, refId=kcid)
            if rmsd < rmsd_cutoff:
                keep = False
                break
        if keep:
            kept.append(cid)
    return kept


def generate_conformers_from_inchi(inchi, rmsd_cutoff=0.15):

    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        raise ValueError(f"Could not parse InChI: {inchi!r}")
    mol = Chem.AddHs(mol)

    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_confs = max(5, n_rot * 5)

    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=42))
    if not conf_ids:
        raise RuntimeError("Embedding produced no conformers")

    conf_ids = prune_conformers_by_rmsd_keep_order(mol, conf_ids, rmsd_cutoff=rmsd_cutoff)


    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    energies = []
    for cid in conf_ids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        if ff is not None:
            ff.Minimize()
            energies.append((cid, ff.CalcEnergy()))

    if energies:
        energies.sort(key=lambda x: x[1])
        ordered_conf_ids = [cid for cid, e in energies]
        energy_dict = {cid: e for cid, e in energies}
    else:
        ordered_conf_ids = conf_ids


    final_conf_ids = prune_conformers_by_rmsd_keep_order(mol, ordered_conf_ids, rmsd_cutoff=rmsd_cutoff)

    return mol, final_conf_ids
