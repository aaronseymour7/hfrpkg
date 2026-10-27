"""How often does reaction generation succeed on a diverse set of molecules?  (generation only; no accuracy claim)

Every neutral, closed-shell species of the bundled ATcT table that contains only C, H, N, O, F, Cl and
has 3-12 heavy atoms is run through the isodesmic, hypohomodesmotic and homodesmotic generators.
Each scheme is run twice, with and without the trivial identity reaction (input -> input), to separate
"no balanced reaction exists" from "only the identity reaction exists" (the molecule is its own smallest
fragment).  Records the outcome, the solve time and, for errors, the exception.  Output: data/feasibility.csv
"""
import csv
import importlib.resources as resources
import signal
import time

from rdkit import Chem, RDLogger

import hfrpkg.core as core
from hfrpkg.core import Homodesmotic, Hypohomodesmotic, Isodesmic

RDLogger.DisableLog("rdApp.*")
ALLOWED = {"C", "H", "N", "O", "F", "Cl"}
FUNCS = (("ID", Isodesmic), ("HY", Hypohomodesmotic), ("HO", Homodesmotic))
GROUPS = {
    "aromatic": "a", "nitrile": "C#N", "amide": "C(=O)[NX3]", "nitro": "[N+](=O)[O-]",
    "N-containing": "[#7]", "halogen": "[F,Cl]", "carbonyl": "[CX3]=O", "ring": "[R]",
}


class Timeout(Exception):
    pass


def _alarm(*_):
    raise Timeout()


def main(limit_s=60):
    signal.signal(signal.SIGALRM, _alarm)
    text = resources.files("hfrpkg.data").joinpath("ATcT_lib.txt").read_text(encoding="utf-8")
    pats = {k: Chem.MolFromSmarts(v) for k, v in GROUPS.items()}
    out = []
    for line in text.splitlines():
        p = line.split("\t")
        if len(p) < 7:
            continue
        m = Chem.MolFromSmiles(p[4]) or Chem.MolFromInchi(p[3])
        if m is None or Chem.GetFormalCharge(m) != 0:
            continue
        mh = Chem.AddHs(m)
        if any(a.GetNumRadicalElectrons() for a in mh.GetAtoms()) or {a.GetSymbol() for a in mh.GetAtoms()} - ALLOWED:
            continue
        if not 3 <= m.GetNumHeavyAtoms() <= 12:
            continue
        flags = [k for k, q in pats.items() if m.HasSubstructMatch(q)]
        for scheme, fn in FUNCS:
            t = time.time()
            outcome, detail = "", ""
            for exclude in (True, False):
                core.EXCLUDE_IDENTITY = exclude
                signal.alarm(limit_s)
                try:
                    _, _, status = fn(Chem.AddHs(m), None, None, None, None)
                except Timeout:
                    outcome, detail = "timeout", f">{limit_s}s"
                except Exception as exc:  # recorded, not hidden
                    outcome, detail = "error", f"{type(exc).__name__}: {str(exc)[:60]}"
                finally:
                    signal.alarm(0)
                    core.EXCLUDE_IDENTITY = True
                if outcome in ("timeout", "error"):
                    break
                if exclude and status == "Optimal":
                    outcome = "reaction"
                    break
                if not exclude:
                    outcome = "identity only" if status == "Optimal" else "no reaction"
            out.append([p[3], Chem.MolToSmiles(m), m.GetNumHeavyAtoms(), ";".join(flags), scheme, outcome, f"{time.time() - t:.3f}", detail])
    with open("data/feasibility.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["inchi", "smiles", "heavy_atoms", "groups", "scheme", "outcome", "seconds", "detail"])
        w.writerows(out)
    print(len(out) // 3, "molecules scanned")


if __name__ == "__main__":
    main()
