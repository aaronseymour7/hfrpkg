import argparse
import glob
import os
import sys

from hfrpkg.common import read_index, resolve_software, write_index
from hfrpkg.conformers import select_lowest
from hfrpkg.run_single import add_job_args, run_jobs


def make_spec(method, basis, extension, folder=".", submit=True, **job_kwargs):
    """Write single-point inputs into <folder>/spec from the optimised logs, then submit."""
    from AaronTools.geometry import Geometry
    from AaronTools.theory import Theory
    from AaronTools.theory.job_types import SinglePointJob

    index = read_index(folder)
    _, opt_out = index.extensions
    sp_software, sp_in, _ = resolve_software(extension)

    if not glob.glob(os.path.join(folder, "*" + opt_out)):
        print(f"No {opt_out} files found in {folder}.")
        sys.exit(1)

    spec_dir = os.path.join(folder, "spec")
    os.makedirs(spec_dir, exist_ok=True)
    write_index(spec_dir, index.level, sp_software, index.input_smiles, index.entries)

    theory = Theory(method=method, basis=basis, job_type=SinglePointJob())
    for name, _, _ in index.entries:
        out = os.path.join(spec_dir, name + sp_in)
        try:
            log, _, info = select_lowest(folder, name, index.conformers.get(name, 1), opt_out)
            Geometry(log).write(outfile=out, theory=theory)
            print(f"Wrote: {out}" + (f"  (from {info['selected']}, lowest of {info['n']} conformers)" if info["n"] > 1 else ""))
        except Exception as exc:
            print(f"Error processing {name}: {exc}")

    if submit:
        run_jobs(spec_dir, **job_kwargs)


def main_cli():
    parser = argparse.ArgumentParser(description="Generate and submit single-point inputs from optimised logs")
    parser.add_argument("folder", nargs="?", default=".")
    parser.add_argument("--m", "--method", dest="method", default="m06-2x", help="default: m06-2x")
    parser.add_argument("--b", "--basis", dest="basis", default="def2tzvp", help="default: def2tzvp")
    parser.add_argument("--s", choices=["g", "o", "p"], default="g", help="g/o/p (default g)")
    add_job_args(parser)
    a = parser.parse_args()
    make_spec(a.method, a.basis, a.s, a.folder, walltime=a.walltime, processors=a.cores,
              memory=a.mem, rerun=a.rerun)


if __name__ == "__main__":
    main_cli()
