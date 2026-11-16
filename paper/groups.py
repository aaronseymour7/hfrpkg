"""The three molecule groups used in the manuscript; every benchmark molecule is in exactly one.

    aromatic  contains a benzene ring (including benzo-fused small rings such as indane)
    strained  otherwise contains a 3-5-membered ring, or is a cyclic allene
              (1,2,4-cyclohexatriene, cycloheptatetraene)
    simple    acyclic molecules and unstrained 6-7-membered rings

`cls` is the class column of data/benchmark_errors.csv (unstrained / strained / ring).
"""
from rdkit import Chem

GROUPS = ("simple", "strained", "aromatic")
BENZENE = Chem.MolFromSmarts("c1ccccc1")
CYCLIC_ALLENES = {"C1=CC=CCC=1", "C1=CC=CC=CC=1"}


def group(smiles, cls):
    m = Chem.MolFromSmiles(smiles)
    if m.HasSubstructMatch(BENZENE):
        return "aromatic"
    if cls == "strained" or smiles in CYCLIC_ALLENES:
        return "strained"
    return "simple"


def n_benzene(smiles):
    return len(Chem.MolFromSmiles(smiles).GetSubstructMatches(BENZENE))
