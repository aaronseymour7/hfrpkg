"""Collect the benchmark workbooks into one tidy CSV (data/benchmark_errors.csv).

Usage: python export_data.py OPT_XLSX SP_XLSX
Errors are in kcal/mol relative to ATcT (0 K). The last 18 rows of each sheet are the
strained (3-5-membered ring) species; the first 54 are acyclic.
"""
import csv
import sys

import openpyxl
from rdkit import Chem

LEVELS = {  # sheet -> (stage, level)
    "B3LYP631GD": ("opt", "B3LYP/6-31+G(d)"),
    "M062X631GD": ("opt", "M06-2X/6-31+G(d)"),
    "wB97DEF2tzvp": ("opt", "wB97X-D/def2-TZVP"),
    "SPB3L": ("sp", "B3LYP/6-31+G(d)"),
    "SPM06": ("sp", "M06-2X/6-31+G(d)"),
    "SPwB97": ("sp", "wB97X-D/def2-TZVP"),
}


def main(opt_xlsx, sp_xlsx, out="data/benchmark_errors.csv"):
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["inchi", "smiles", "class", "atct_kcal", "stage", "geometry_level", "ID", "HY", "HO"])
        for path in (opt_xlsx, sp_xlsx):
            wb = openpyxl.load_workbook(path, data_only=True)
            for ws in wb:
                if ws.title not in LEVELS:
                    continue
                stage, level = LEVELS[ws.title]
                rows = [r[:5] for r in ws.iter_rows(min_row=2, values_only=True) if r[0]]
                for i, (inchi, atct, idv, hy, ho) in enumerate(rows):
                    smi = Chem.MolToSmiles(Chem.MolFromInchi(inchi))
                    w.writerow([inchi, smi, "strained" if i >= 54 else "unstrained", repr(float(atct)), stage, level,
                                repr(float(idv)), repr(float(hy)), repr(float(ho))])


if __name__ == "__main__":
    main(*sys.argv[1:3])
