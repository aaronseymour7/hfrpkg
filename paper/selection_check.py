"""Which survey species satisfy the benchmark's selection criteria?

(1) Decomposable: a non-trivial reaction exists under each of the three schemes (isodesmic,
    hypohomodesmotic, homodesmotic), i.e. the molecule is not its own smallest fragment.
(2) Verifiable: every reference species in those reactions has an ATcT enthalpy of formation, so that the
    computed value of the target can be checked.

The reactions are generated with the solver formulation used for the published benchmark (no objective,
core.USE_OBJECTIVE = False; identity not excluded, core.EXCLUDE_IDENTITY = False), which paper/objective_check.py
shows to agree with the current formulation for all 72 benchmark molecules except four hypohomodesmotic cases.
Every qualifying species is listed with its status: in the benchmark, or excluded, with the reason.
Species list: data/feasibility.csv.  Output: data/selection_check.csv
"""
import csv

from rdkit import Chem, RDLogger

import hfrpkg.core as core
from hfrpkg.common import get_Hf

RDLogger.DisableLog("rdApp.*")
SCHEMES = ("Isodesmic", "Hypohomodesmotic", "Homodesmotic")
EXCLUDED = {  # qualifying species left out of the study, and why
    "InChI=1S/C2H4O2/c3-1-2-4/h1-4H/": "E/Z isomer not specified by ATcT",
    "InChI=1S/C3H6O2/c1-3(5)2-4/h2,4-5H,1H3": "E/Z isomer not specified by ATcT",
    "InChI=1S/C2H2Cl2/c3-1-2-4/h1-2H": "stereo-unspecified duplicate of the (Z) isomer in the benchmark",
    "InChI=1S/C2H2F2/c3-1-2-4/h1-2H": "stereo-unspecified duplicate of the (Z) isomer in the benchmark",
    "InChI=1S/2CH4O/c2*1-2/h2*2H,1H3": "dimer (two molecules), not a single species",
    "InChI=1S/C5H4/c1-2-4-3(1)5(1,2)4/h1-4H": "pyramidane: no complete conformer set at wB97X-D/def2-TZVP",
}


def main():
    core.USE_OBJECTIVE, core.EXCLUDE_IDENTITY = False, False
    species = list(dict.fromkeys((r["inchi"], r["smiles"]) for r in csv.DictReader(open("data/feasibility.csv", encoding="utf-8"))))
    bench = {r["inchi"] for r in csv.DictReader(open("data/benchmark_errors.csv", encoding="utf-8"))}
    out = []
    for inchi, smiles in species:
        decomposable, missing = True, set()
        for name in SCHEMES:
            m = Chem.AddHs(Chem.MolFromSmiles(smiles))
            result = getattr(core, name)(m, None, None, None, None)
            if len(result) < 3 or result[2] != "Optimal":
                decomposable = False
                continue
            rhs, lhs, _ = result
            target = Chem.MolToInchi(m)
            if any(Chem.MolToInchi(x) == target for x, _c in rhs):
                decomposable = False
                continue
            missing |= {Chem.MolToSmiles(Chem.RemoveHs(x)) for x, _c in list(lhs) + list(rhs)
                        if Chem.MolToInchi(x) != target and get_Hf(Chem.MolToInchi(x)) is None}
        qualifies = decomposable and not missing
        status = "benchmark" if inchi in bench else "excluded: " + EXCLUDED[inchi] if qualifies else ""
        out.append([inchi, smiles, inchi in bench, decomposable, qualifies, ";".join(sorted(missing)), status])
    core.USE_OBJECTIVE, core.EXCLUDE_IDENTITY = True, True
    with open("data/selection_check.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["inchi", "smiles", "in_benchmark", "decomposable_all3", "decomposable_and_all_references_in_ATcT",
                    "references_without_ATcT", "status"])
        w.writerows(out)
    n, dec, ok = len(out), sum(r[3] for r in out), sum(r[4] for r in out)
    print(f"{n} survey species; (1) decomposable under all three: {dec}; (1)+(2) also every reference has ATcT: {ok}")
    print(f"benchmark molecules satisfying both: {sum(r[2] and r[4] for r in out)} of {sum(r[2] for r in out)}; "
          f"qualifying species not in the benchmark (excluded): {sum((not r[2]) and r[4] for r in out)}")


if __name__ == "__main__":
    main()
