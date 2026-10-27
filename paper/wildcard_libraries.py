"""The wildcard fragment libraries of each level, extracted from hfrpkg/core.py (data/wildcard_libraries.csv).

Each pattern is a generic fragment in which every [*] stands for any atom.  The generator maps a pattern
onto the (kekulized) target, replaces each [*] by the element of the matched atom, adds hydrogens, and
keeps the distinct species it obtains: "products" are candidate right-hand-side species, "reactants" the
candidate left-hand-side co-reactants.  Level numbers are those used in core.get_fragments.

Level 2 is the isodesmic scheme.  Isodesmic does not map this library: core.Isodesmic builds the same
candidates directly, one species X-X', X=X' or X#X' (hydrogens added) for each atom-pair/bond-order type
present in the target, the hydride of each element present, and H2.  The patterns below are those generic
forms; the H2 co-reactant is added explicitly by core.Isodesmic and is therefore appended here.
"""
import csv
import inspect
import re

import hfrpkg.core as core

SCHEME = {1: "isogyric", 2: "isodesmic", 3: "hypohomodesmotic", 4: "homodesmotic"}


def main():
    src = inspect.getsource(core.get_fragments)
    rows, level, role = [], None, None
    for line in src.splitlines():
        m = re.match(r"\s*if level\s*==\s*(\d+)\s*:", line)
        if m:
            level = int(m.group(1))
        m = re.match(r"\s*wildcard_(reactants|products)\s*=\s*\[", line)
        if m:
            role = m.group(1)
        for kind, pat in re.findall(r'Chem\.MolFrom(Smiles|Smarts)\("([^"]+)"\)', line):
            if level and role:
                rows.append([level, SCHEME[level], role, pat, kind.upper(), pat.count("*")])
    rows.insert(next(i for i, r in enumerate(rows) if r[0] == 2 and r[2] == "products"), [2, SCHEME[2], "reactants", "[H][H]", "SMILES", 0])
    with open("data/wildcard_libraries.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["level", "scheme", "role", "pattern", "syntax", "wildcard_atoms"])
        w.writerows(rows)
    for lv in sorted(SCHEME):
        n = [r for r in rows if r[0] == lv]
        print(f"level {lv} ({SCHEME[lv]}): {sum(r[2] == 'reactants' for r in n)} reactant and {sum(r[2] == 'products' for r in n)} product patterns")


if __name__ == "__main__":
    main()
