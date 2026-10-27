import csv
import glob
import os

from hfrpkg.multi.fill import fill_logs
from hfrpkg.multi.grab_unique import numeric_key
from hfrpkg.spec_compute import spec_compute
from hfrpkg.write_data import CSV_HEADER, write_summary


def main():
    fill_logs("unique_files", subdir="spec")
    rows = []
    with open("sp_reaction_summaries.txt", "w", encoding="utf-8") as rxn_fout:
        for mhfr in sorted(glob.glob("*.mhfr"), key=numeric_key):
            data = spec_compute(mhfr)
            if data is None:
                continue
            path = os.path.join(mhfr, "spec", "sp_summary.txt")
            rows.append(write_summary(data, os.path.join(mhfr, "spec"), "sp_summary.txt"))
            with open(path, encoding="utf-8") as sf:
                rxn_fout.write(f"=== {path} ===\n{sf.read()}\n")

    with open("sp_enthalpies_summary.csv", "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(CSV_HEADER)
        writer.writerows(rows)


def main_cli():
    main()


if __name__ == "__main__":
    main_cli()
