"""Does the ILP objective change the reactions?  Compares three formulations on the benchmark molecules.

    python objective_check.py path/to/legacy/hfrpkg/core.py

legacy core.py : the benchmark-era code (no objective).  It is run as-is and with the stray
                 `CC=[NH]` hypohomodesmotic option removed; both are compared with the current
                 code (minimum-species objective, identity reaction excluded).
Writes data/objective_check.csv.
"""
import csv
import importlib.util
import re
import sys
import tempfile

from rdkit import Chem, RDLogger

import hfrpkg.core as NEW

RDLogger.DisableLog("rdApp.*")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main(legacy_path):
    src = open(legacy_path, encoding="utf-8").read()
    nohack = re.sub(r"[ \t]*rhs_options\.append\(Chem\.MolFromSmiles\('CC=\[NH\]'\)\)\n", "", src)
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(nohack)
    variants = {"legacy": load(legacy_path, "legacy"), "legacy_no_imine": load(fh.name, "legacy2"), "current": NEW}
    mols = []
    for r in csv.DictReader(open("data/benchmark_errors.csv", encoding="utf-8")):
        if r["smiles"] not in mols:
            mols.append(r["smiles"])
    key = lambda g: sorted((c, Chem.MolToSmiles(Chem.RemoveHs(m))) for m, c in g)
    fmt = lambda x: " + ".join((f"{c} " if c > 1 else "") + s for c, s in x[0]) + " -> " + " + ".join((f"{c} " if c > 1 else "") + s for c, s in x[1])
    summary = {}
    with open("data/objective_check.csv", "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["smiles", "scheme"] + list(variants))
        for smi in mols:
            for fn, k in (("Isodesmic", "ID"), ("Hypohomodesmotic", "HY"), ("Homodesmotic", "HO")):
                res = {}
                for tag, mod in variants.items():
                    rhs, lhs, st = getattr(mod, fn)(Chem.AddHs(Chem.MolFromSmiles(smi)), None, None, None, None)
                    res[tag] = (key(lhs), key(rhs)) if st == "Optimal" else None
                w.writerow([smi, k] + [fmt(res[t]) if res[t] else "infeasible" for t in variants])
                summary.setdefault(k, [0, 0, 0])
                summary[k][0] += res["legacy_no_imine"] == res["current"]
                summary[k][1] += res["legacy"] == res["current"]
                summary[k][2] += 1
    for k, (a, b, n) in summary.items():
        print(f"{k}: legacy-without-imine == current for {a}/{n}; legacy == current for {b}/{n}")


if __name__ == "__main__":
    main(sys.argv[1])
