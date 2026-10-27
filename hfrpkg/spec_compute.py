#!/usr/bin/env python3
"""Single-point post-processing for one folder (needs <folder>/spec/*.log)."""

import argparse
import os

from hfrpkg.compute_folder import compute_spec_folder
from hfrpkg.write_data import write_summary


def spec_compute(folder="."):
    data = compute_spec_folder(folder)
    if data is None:
        return None
    write_summary(data, os.path.join(folder, "spec"), "sp_summary.txt")
    return data


def main_cli():
    parser = argparse.ArgumentParser(description="Compute dHf from single-point (spec/) results")
    parser.add_argument("--f", "--folder", dest="folder", default=".", help="reaction folder (default: .)")
    args = parser.parse_args()
    spec_compute(args.folder)


if __name__ == "__main__":
    main_cli()
