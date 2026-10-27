from collections import Counter

import pytest
from rdkit import Chem

from hfrpkg.core import Homodesmotic, Hypohomodesmotic, Isodesmic, Isogyric


def atoms(group):
    total = Counter()
    for mol, coeff in group:
        for a in Chem.AddHs(mol).GetAtoms():
            total[a.GetSymbol()] += coeff
    return total


@pytest.mark.parametrize("fn", [Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic])
@pytest.mark.parametrize("smiles", ["CCCCO", "c1ccccc1", "CC(C)(C)CC"])
def test_reaction_conserves_atoms_and_is_not_trivial(fn, smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rhs, lhs, status = fn(mol, None, None, None, None)
    assert status == "Optimal"
    assert atoms(lhs) == atoms(rhs)
    input_inchi = Chem.MolToInchi(mol)
    # the input appears exactly once on the left with coefficient 1 ...
    assert [c for m, c in lhs if Chem.MolToInchi(m) == input_inchi] == [1]
    # ... and never as its own product (identity reaction)
    assert all(Chem.MolToInchi(m) != input_inchi for m, _ in rhs)


def test_infeasible_is_reported_not_raised():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    rhs, lhs, status = Homodesmotic(mol, None, None, None, None)
    assert status != "Optimal" and rhs == [] and lhs == []


def test_identity_switch_separates_no_reaction_from_identity_only():
    """Ethanol is its own smallest homodesmotic fragment: only the trivial reaction exists."""
    import hfrpkg.core as core

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert Homodesmotic(mol, None, None, None, None)[2] != "Optimal"  # trivial reaction excluded
    core.EXCLUDE_IDENTITY = False
    try:
        rhs, lhs, status = Homodesmotic(mol, None, None, None, None)
    finally:
        core.EXCLUDE_IDENTITY = True
    assert status == "Optimal"
    inchi = Chem.MolToInchi(mol)
    assert [Chem.MolToInchi(m) for m, _ in rhs] == [inchi]  # identity: ethanol -> ethanol


def test_solver_switches_reproduce_the_benchmark_formulation():
    """USE_OBJECTIVE=False with EXCLUDE_IDENTITY=False is the formulation that generated the benchmark."""
    import hfrpkg.core as core

    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCO"))
    default = Homodesmotic(mol, None, None, None, None)
    core.USE_OBJECTIVE, core.EXCLUDE_IDENTITY = False, False
    try:
        legacy = Homodesmotic(mol, None, None, None, None)
    finally:
        core.USE_OBJECTIVE, core.EXCLUDE_IDENTITY = True, True
    def key(r):
        side = lambda g: sorted((c, Chem.MolToSmiles(Chem.RemoveHs(m))) for m, c in g)
        return side(r[1]), side(r[0])

    assert default[2] == legacy[2] == "Optimal" and key(default) == key(legacy)  # unique reaction: formulation does not matter
