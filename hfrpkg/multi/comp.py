import csv
import glob
import os

from hfrpkg.compute_folder import compute_folder
from hfrpkg.multi.grab_unique import numeric_key
from hfrpkg.write_data import CSV_HEADER, write_conformer_selection, write_single_reaction


def main():
    rows = []
    with open("Reactions_summary.txt", "w", encoding="utf-8") as rxn_fout:
        for mhfr in sorted(glob.glob("*.mhfr"), key=numeric_key):
            if not os.path.isdir(mhfr):
                continue
            data = compute_folder(mhfr)
            if data is None:
                continue
            rows.append(write_single_reaction(data, mhfr))
            write_conformer_selection(data, mhfr)
            with open(os.path.join(mhfr, "opt_summary.txt"), encoding="utf-8") as sf:
                rxn_fout.write(f"=== {mhfr}/opt_summary.txt ===\n{sf.read()}\n")

    with open("enthalpies_summary.csv", "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(CSV_HEADER)
        writer.writerows(rows)
