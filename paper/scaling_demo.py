"""Does the homodesmotic reference set grow with molecule size?  (It should not.)

Generates the homodesmotic reaction for polyenes, polycyclic aromatics and C60 and records
the number of distinct reference species, their total stoichiometry and the ILP solve time.
Output: data/scaling.csv
"""
import csv
import time
from collections import Counter

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from hfrpkg.common import get_Hf
from hfrpkg.core import Homodesmotic
from hfrpkg.utils import homodesmotic_count

C60 = ("C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4"
       "C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23")
CASES = (
    [("polyene", f"C{2 * n + 2}", "C=C" + "C=C" * n) for n in range(2, 20, 2)]
    + [("PAH", "benzene", "c1ccccc1"), ("PAH", "naphthalene", "c1ccc2ccccc2c1"),
       ("PAH", "anthracene", "c1ccc2cc3ccccc3cc2c1"), ("PAH", "pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
       ("PAH", "coronene", "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"), ("fullerene", "C60", C60)]
)


def conserved(lhs, rhs):
    """Independent check: bond hybridization types and hydrogen environments balance."""
    def side(group):
        b, h = Counter(), Counter()
        for m, c in group:
            bb, hh = homodesmotic_count(m)
            for k, v in bb.items():
                b[k] += c * v
            for k, v in hh.items():
                h[k] += c * v
        return b, h
    return side(lhs) == side(rhs)


def smi(m):
    return Chem.MolToSmiles(Chem.RemoveHs(m))


def main(out="data/scaling.csv"):
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["family", "name", "formula", "heavy_atoms", "solve_s", "distinct_refs", "total_ref_molecules",
                    "all_refs_in_ATcT", "conserved", "reaction"])
        for fam, name, s in CASES:
            m = Chem.MolFromSmiles(s)
            inp = Chem.MolToInchi(Chem.AddHs(m))
            t = time.time()
            rhs, lhs, status = Homodesmotic(Chem.AddHs(m), None, None, None, None)
            dt = time.time() - t
            if status != "Optimal":
                print(name, status)
                continue
            refs = [(c, x) for g in (lhs, rhs) for x, c in [(a, b) for a, b in g] if Chem.MolToInchi(x) != inp]
            atct = all(get_Hf(Chem.MolToInchi(x)) is not None for _, x in refs)
            fmt = lambda g: " + ".join((f"{c} " if c > 1 else "") + (name if Chem.MolToInchi(x) == inp else smi(x)) for x, c in g)
            w.writerow([fam, name, CalcMolFormula(m), m.GetNumAtoms(), f"{dt:.3f}", len(refs), sum(c for c, _ in refs),
                        atct, conserved(lhs, rhs), f"{fmt(lhs)} -> {fmt(rhs)}"])
            print(f"{name:12s} atoms={m.GetNumAtoms():3d} refs={len(refs)} total={sum(c for c, _ in refs):4d} {dt:.2f}s conserved={conserved(lhs, rhs)}")


if __name__ == "__main__":
    main()
