"""ATcT uncertainties of the benchmark molecules and of the reference species used in their reactions.

Writes data/atct_uncertainty.csv (kcal/mol).  Reference uncertainties are propagated in quadrature,
sigma_ref = sqrt(sum_j (nu_j * sigma_j)^2) / nu_M, which ignores the covariance between ATcT values
(ATcT values share data, so this is a conservative-to-approximate estimate, not a rigorous bound).
The reactions are the ones used in the benchmark: current code, except for the hypohomodesmotic
reactions that depend on the ILP objective, which use the benchmark (no-objective) solution.
"""
import csv
import importlib.resources as resources
import re

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
KJ_PER_KCAL = 4.184


def library():
    text = resources.files("hfrpkg.data").joinpath("ATcT_lib.txt").read_text(encoding="utf-8")
    out = {}
    for line in text.splitlines():
        p = line.split("\t")
        if len(p) >= 7:
            out[p[3]] = (float(p[5]), float(p[6]))  # kJ/mol: dHf(0 K), uncertainty
    return out


def parse(reaction):
    lhs, rhs = reaction.split(" -> ")
    def side(s):
        terms = []
        for t in s.split(" + "):
            m = re.match(r"^(\d+) (.+)$", t)
            terms.append((int(m.group(1)), m.group(2)) if m else (1, t))
        return terms
    return side(lhs), side(rhs)


def inchi(smiles):
    return Chem.MolToInchi(Chem.AddHs(Chem.MolFromSmiles(smiles)))


def main():
    lib = library()
    rows = list(csv.DictReader(open("data/objective_check.csv", encoding="utf-8")))
    targets = {}
    for r in csv.DictReader(open("data/benchmark_errors.csv", encoding="utf-8")):
        targets.setdefault(r["inchi"], (r["smiles"], float(r["atct_kcal"])))
    sig_ref = {}
    for r in rows:
        col = "legacy_no_imine" if r["scheme"] == "HY" else "current"
        lhs, rhs = parse(r[col])
        target = inchi(r["smiles"])
        var = 0.0
        for nu, smi in lhs + rhs:
            i = inchi(smi)
            if i == target:
                continue
            var += (nu * lib[i][1]) ** 2
        sig_ref[(target, r["scheme"])] = var ** 0.5 / KJ_PER_KCAL
    with open("data/atct_uncertainty.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["inchi", "smiles", "atct_kcal", "sigma_target_kcal", "sigma_refs_ID", "sigma_refs_HY", "sigma_refs_HO"])
        for i, (smi, h) in targets.items():
            w.writerow([i, smi, f"{h:.4f}", f"{lib[i][1] / KJ_PER_KCAL:.4f}"] + [f"{sig_ref[(i, k)]:.4f}" for k in ("ID", "HY", "HO")])


if __name__ == "__main__":
    main()
