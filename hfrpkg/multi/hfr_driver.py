#!/usr/bin/env python3
import argparse

reaction_map = {
    "1": "isogyric",
    "2": "isodesmic",
    "3": "hypohomodesmotic",
    "4": "homodesmotic",
}


def main(input_file, method=None, basis=None, software="g", conformers=None):
    from hfrpkg.runner import run_reaction

    with open(input_file, "r") as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                print(f"[SKIP] Line {i} malformed: {line}")
                continue
            level, smiles = parts
            reaction_type = reaction_map.get(level)
            if not reaction_type:
                print(f"[SKIP] Invalid level {level} on line {i}")
                continue
            try:
                run_reaction(
                    action_type="write",
                    reaction_type=reaction_type,
                    input_smiles=smiles,
                    outfolder=f"{i}.mhfr",
                    method=method,
                    basis=basis,
                    extension=software,
                    conformers=conformers,
                )
            except Exception as exc:
                print(f"[ERROR] Line {i} failed: {exc}")


def main_cli():
    parser = argparse.ArgumentParser(description="Write reaction folders for a list of molecules")
    parser.add_argument("input_file", help="text file: '<level 1-4> <SMILES>' per line")
    parser.add_argument("--m", help="METHOD")
    parser.add_argument("--b", help="BASIS")
    parser.add_argument("--s", choices=["g", "o", "p"], default="g", help="g/o/p (default g)")
    from hfrpkg.cli import add_conformer_args, conformer_settings

    add_conformer_args(parser)
    a = parser.parse_args()

    main(a.input_file, method=a.m, basis=a.b, software=a.s, conformers=conformer_settings(a))


if __name__ == "__main__":
    main_cli()
