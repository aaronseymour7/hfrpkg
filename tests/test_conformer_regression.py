"""Conformer sets for the 72 benchmark molecules under the RDKit used for the published data.

Embedding depends on the RDKit version, so this only runs under RDKit 2024.03.3 (skipped otherwise;
CI has a job that installs it).  The reference counts and MMFF94s energies were produced by this
package and verified identical to the original makeConf.py under that version.
"""
import csv
import os

import pytest
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from hfrpkg.conformers import generate_conformers

DATA = os.path.join(os.path.dirname(__file__), "..", "paper", "data", "conformer_counts_rdkit2024.03.3.csv")

pytestmark = pytest.mark.skipif(rdkit.__version__ != "2024.03.3", reason="needs RDKit 2024.03.3 (the benchmark version)")


def test_benchmark_conformer_sets_are_reproduced():
    rows = list(csv.DictReader(open(DATA, encoding="utf-8")))
    assert len(rows) == 72
    for r in rows:
        mol, ids = generate_conformers(r["inchi"])
        assert len(ids) == int(r["n_conformers"]), r["smiles"]
        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        energies = [AllChem.MMFFGetMoleculeForceField(mol, mp, confId=c).CalcEnergy() for c in ids]
        assert energies == pytest.approx([float(x) for x in r["mmff94s_energies_kcal_per_mol"].split(";")], abs=2e-3), r["smiles"]
