#!/usr/bin/env python3

import argparse


def parse_smiles_list(arg):
    from rdkit import Chem

    if not arg:
        return None
    raw = ",".join(arg) if isinstance(arg, list) else arg
    mols = []
    for smi in raw.split(","):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise argparse.ArgumentTypeError(f"Invalid SMILES: {smi}")
        mols.append(mol)
    return mols


def build_parser():
    parser = argparse.ArgumentParser(description="Run reaction balancing (HFR)")
    parser.add_argument("action_type", choices=["write", "view", "count"], help="Action to perform")
    parser.add_argument(
        "reaction_type",
        choices=["isogyric", "isodesmic", "homodesmotic", "hypohomodesmotic"],
        help="Reaction type",
    )
    parser.add_argument("input", type=str, help="Input molecule SMILES")
    parser.add_argument("--lhs", nargs="+", help="Required LHS molecules (comma-separated or list)")
    parser.add_argument("--rhs", nargs="+", help="Required RHS molecules (comma-separated or list)")
    parser.add_argument("--substruct", nargs="+", help="Substructures to replace")
    parser.add_argument("--replacement", nargs="+", help="Replacement structures")
    parser.add_argument("--outfolder", help="Output folder (required for write)")
    parser.add_argument("--m", help="Method (default B3LYP)")
    parser.add_argument("--b", help="Basis set (default 6-31G)")
    parser.add_argument(
        "--s", choices=["g", "o", "p"], default="g",
        help="Software: g=Gaussian, o=ORCA, p=Psi4 (default g)",
    )
    add_conformer_args(parser)
    return parser


def add_conformer_args(parser):
    parser.add_argument("--no-conformers", action="store_true", help="embed a single conformer per species")
    parser.add_argument("--min-conformers", type=int, default=5, help="minimum conformers embedded per species (default 5)")
    parser.add_argument("--per-rotor", type=int, default=5, help="conformers per rotatable bond, n = max(min, per_rotor*n_rot) (default 5)")
    parser.add_argument("--conf-rmsd", type=float, default=0.15, help="RMSD pruning cutoff in Angstrom (default 0.15)")
    parser.add_argument("--conf-seed", type=int, default=42, help="embedding random seed (default 42)")


def conformer_settings(args):
    from hfrpkg.conformers import ConformerSettings

    return ConformerSettings(enabled=not args.no_conformers, min_conformers=args.min_conformers,
                             per_rotor=args.per_rotor, rmsd_cutoff=args.conf_rmsd, seed=args.conf_seed)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.action_type == "write" and not args.outfolder:
        parser.error("--outfolder is required for the write action")

    from hfrpkg.runner import run_reaction  # heavy imports only after argparse

    run_reaction(
        action_type=args.action_type,
        reaction_type=args.reaction_type,
        input_smiles=args.input,
        lhs=parse_smiles_list(args.lhs),
        rhs=parse_smiles_list(args.rhs),
        substruct=parse_smiles_list(args.substruct),
        replacement=parse_smiles_list(args.replacement),
        outfolder=args.outfolder,
        method=args.m,
        basis=args.b,
        extension=args.s,
        conformers=conformer_settings(args),
    )


if __name__ == "__main__":
    main()
