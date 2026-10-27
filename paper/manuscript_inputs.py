"""Derive the inputs of build_manuscript.js from the data in this repository.

    python manuscript_inputs.py   ->  data/stats.json, data/mols.json, data/hy_ambiguous.json

stats.json         error statistics per stage / subset / level / scheme (Tables 1, 2, S1)
mols.json          the 77 benchmark molecules with class, formula, ring sizes (Table S2)
hy_ambiguous.json  hypohomodesmotic reactions that depend on the ILP objective (Table S4)
wildcards.json     wildcard fragment libraries (Table S5)
feasibility_summary.json  applicability survey (feasibility.csv)
"""
import csv
import json
import statistics

from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
KEY = {  # (stage, geometry level) -> sheet key used by the manuscript builder
    ("opt", "B3LYP/6-31+G(d)"): "B3LYP631GD", ("opt", "M06-2X/6-31+G(d)"): "M062X631GD", ("opt", "wB97X-D/def2-TZVP"): "wB97DEF2tzvp",
    ("sp", "B3LYP/6-31+G(d)"): "SPB3L", ("sp", "M06-2X/6-31+G(d)"): "SPM06", ("sp", "wB97X-D/def2-TZVP"): "SPwB97",
}
SCHEMES = ("ID", "HY", "HO")


def summarize(v):
    n = len(v)
    me = sum(v) / n
    return dict(n=n, rmse=(sum(x * x for x in v) / n) ** 0.5, mae=sum(abs(x) for x in v) / n, me=me, std=statistics.stdev(v))


def main():
    rows = list(csv.DictReader(open("data/benchmark_errors.csv", encoding="utf-8")))
    stats = {}
    for stage in ("opt", "sp"):
        for subset in ("unstrained", "strained", "combined"):
            sel = [r for r in rows if r["stage"] == stage and (subset == "combined" or r["class"] == subset)]
            for k in SCHEMES:
                stats[f"{stage}|{subset}|pooled|{k}"] = summarize([float(r[k]) for r in sel])
                for (st, lv), key in KEY.items():
                    if st == stage:
                        stats[f"{stage}|{subset}|{key}|{k}"] = summarize([float(r[k]) for r in sel if r["geometry_level"] == lv])
    json.dump(stats, open("data/stats.json", "w"), indent=1)

    unc = {r["inchi"]: r for r in csv.DictReader(open("data/atct_uncertainty.csv", encoding="utf-8"))}
    mols, seen = [], set()
    for r in rows:
        if r["inchi"] in seen:
            continue
        seen.add(r["inchi"])
        m = Chem.MolFromInchi(r["inchi"])
        u = unc[r["inchi"]]
        mols.append(dict(i=len(mols) + 1, cls=r["class"], smiles=Chem.MolToSmiles(m), formula=CalcMolFormula(m),
                         atct=float(r["atct_kcal"]), rings=sorted(len(x) for x in m.GetRingInfo().AtomRings()), inchi=r["inchi"],
                         sigma=float(u["sigma_target_kcal"]), sigma_ref=float(u["sigma_refs_HO"])))
    json.dump(mols, open("data/mols.json", "w"), indent=1)

    amb = [[r["smiles"], r["legacy_no_imine"], r["current"]] for r in csv.DictReader(open("data/objective_check.csv"))
           if r["scheme"] == "HY" and r["legacy_no_imine"] != r["current"]]
    json.dump(amb, open("data/hy_ambiguous.json", "w"), indent=1)
    # wildcard libraries (Table S5) and applicability survey (feasibility.csv)
    json.dump([[int(r["level"]), r["scheme"], r["role"], r["pattern"], r["syntax"]] for r in csv.DictReader(open("data/wildcard_libraries.csv"))],
              open("data/wildcards.json", "w"), indent=1)
    feas = list(csv.DictReader(open("data/feasibility.csv", encoding="utf-8")))
    groups = [("All species", None), ("Contain N", "N-containing"), ("Contain F or Cl", "halogen"), ("Contain C=O", "carbonyl"),
              ("Contain a ring", "ring"), ("Nitriles", "nitrile"), ("Amides", "amide"), ("Aromatic", "aromatic")]
    table = []
    for label, g in groups:
        sel = [r for r in feas if g is None or g in r["groups"].split(";")]
        n = len({r["inchi"] for r in sel})
        def cnt(scheme, out):
            return sum(r["scheme"] == scheme and r["outcome"] == out for r in sel)
        table.append([label, n, cnt("ID", "reaction")] + [cnt(sc, o) for sc in ("HY", "HO") for o in ("reaction", "identity only", "no reaction")])
    json.dump(table, open("data/feasibility_summary.json", "w"), indent=1)
    print(f"stats: {len(stats)} entries, molecules: {len(mols)}, objective-dependent HY reactions: {len(amb)}")


if __name__ == "__main__":
    main()
