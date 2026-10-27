"""Conformer generation, and lowest-energy conformer selection.

Ported from the makeConf.py / extractConfs.py scripts used to produce the published benchmark:

1. Start from the species' InChI, add hydrogens.
2. Embed n = max(5, 5 * n_rotatable_bonds) conformers with RDKit's EmbedMultipleConfs (seed 42).
3. Prune by best-fit all-atom RMSD (0.15 A), keeping the earliest conformer of each cluster.
4. Minimize every survivor with MMFF94s and order by force-field energy.
5. Prune again (RMSD 0.15 A) in that order, so the lowest-energy member of each cluster survives.
6. Optimize every conformer with the DFT method of the study; the one with the lowest E + ZPVE
   from its frequency job is used for the thermochemistry and for single-point refinement.

Layout: species left with one conformer keep the legacy file stem ``R1_2``; species with several use
``R1_2-c01 ... R1_2-cNN`` (c01 = lowest force-field energy) and ``index.txt`` records the count and the
settings. The embedded conformers depend on the RDKit version, which is therefore recorded as well.
"""

import os
from dataclasses import dataclass

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors

from hfrpkg.common import conformer_name


@dataclass(frozen=True)
class ConformerSettings:
    enabled: bool = True  # False: embed a single conformer (same pipeline, n = 1)
    min_conformers: int = 5
    per_rotor: int = 5
    rmsd_cutoff: float = 0.15
    seed: int = 42
    mmff_variant: str = "MMFF94s"

    def count(self, n_rotatable):
        if not self.enabled:
            return 1
        return max(self.min_conformers, n_rotatable * self.per_rotor)

    def describe(self):
        if not self.enabled:
            return f"single conformer (rdkit {rdkit.__version__}, seed {self.seed}, {self.mmff_variant})"
        return (f"rdkit {rdkit.__version__} EmbedMultipleConfs seed={self.seed} "
                f"n=max({self.min_conformers},{self.per_rotor}*n_rot) rmsd={self.rmsd_cutoff} {self.mmff_variant}")


def prune_conformers_by_rmsd_keep_order(mol, candidate_conf_ids, rmsd_cutoff=0.15):
    """Keep a conformer only if it is at least `rmsd_cutoff` (best-fit RMSD, all atoms) from every kept one.

    Note rdMolAlign.GetBestRMS leaves the probe conformer rigidly aligned onto the reference, as in
    the original script; this does not change any energy.
    """
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


def generate_conformers(inchi, settings=None):
    """Return (mol_with_Hs, [conformer ids]) for the species with this InChI, lowest MMFF energy first."""
    settings = settings or ConformerSettings()
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        raise ValueError(f"Could not parse InChI: {inchi!r}")
    mol = Chem.AddHs(mol)

    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_confs = settings.count(n_rot)

    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=settings.seed))
    if not conf_ids:  # the original script stops here; fall back to random starting coordinates
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=settings.seed, useRandomCoords=True))
    if not conf_ids:
        raise RuntimeError(f"Embedding produced no conformers for {inchi}")

    conf_ids = prune_conformers_by_rmsd_keep_order(mol, conf_ids, rmsd_cutoff=settings.rmsd_cutoff)

    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=settings.mmff_variant)
    energies = []
    if mp is not None:  # no MMFF parameters for this species: keep the embedding order
        for cid in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
            if ff is not None:
                ff.Minimize()
                energies.append((cid, ff.CalcEnergy()))

    if energies:
        energies.sort(key=lambda x: x[1])
        ordered = [cid for cid, _ in energies]
    else:
        ordered = conf_ids

    return mol, prune_conformers_by_rmsd_keep_order(mol, ordered, rmsd_cutoff=settings.rmsd_cutoff)


def conformer_logs(directory, name, n, out_ext):
    """Output files of every conformer of `name` (just [name + out_ext] for single-structure species)."""
    if n <= 1:
        return [os.path.join(directory, name + out_ext)]
    return [os.path.join(directory, conformer_name(name, k) + out_ext) for k in range(1, n + 1)]


def select_lowest(directory, name, n, out_ext):
    """Pick the lowest-energy conformer of a species from its finished frequency jobs.

    Returns (log_path, (E_plus_zpve, zpve), info) where info lists the energy of every conformer.
    Raises ValueError if any conformer has no usable energy (as extractConfs.py does for an unfinished
    job): choosing among a partial set could silently return a conformer that is not the lowest.
    """
    from hfrpkg.compute_folder import read_opt_energy  # lazy: avoids a circular import

    logs = conformer_logs(directory, name, n, out_ext)
    results = [read_opt_energy(p) if os.path.exists(p) else (None, None) for p in logs]
    missing = [os.path.basename(p) for p, r in zip(logs, results) if r[0] is None]
    if missing:
        raise ValueError(f"no energy for conformer(s) of {name}: {', '.join(missing)}")
    best = min(range(len(logs)), key=lambda i: results[i][0])
    info = [(os.path.basename(p), r[0]) for p, r in zip(logs, results)]
    return logs[best], results[best], {"selected": os.path.basename(logs[best]), "n": len(logs), "energies": info}
