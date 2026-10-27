#!/usr/bin/env python3
import argparse
import csv

from hfrpkg.compute_folder import compute_folder
from hfrpkg.write_data import CSV_HEADER, write_conformer_selection, write_single_reaction


def main(folder="."):
    data = compute_folder(folder)
    if data is None:
        print("Nothing computed; see messages above.")
        return None
    row = write_single_reaction(data, folder)
    write_conformer_selection(data, folder)
    with open("enthalpy_summary.csv", "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(CSV_HEADER)
        writer.writerow(row)
    return row


def main_cli():
    parser = argparse.ArgumentParser(description="Compute dHf from the finished jobs in a folder")
    parser.add_argument("folder", nargs="?", default=".", help="reaction folder (default: .)")
    main(parser.parse_args().folder)


if __name__ == "__main__":
    main_cli()
