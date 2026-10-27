import argparse

from hfrpkg.multi.compute_multiple import combined_postprocess as compute_multiple
from hfrpkg.multi.unique_spec import make_spec
from hfrpkg.run_single import add_job_args


def multi_spprep(method, basis, extension, **job_kwargs):
    compute_multiple()  # collect optimised logs and summarise the optimisation level first
    make_spec(method, basis, extension, **job_kwargs)


def main_cli():
    parser = argparse.ArgumentParser(description="Generate and submit single-point jobs for all unique molecules")
    parser.add_argument("--m", "--method", dest="method", default="m06-2x", help="default: m06-2x")
    parser.add_argument("--b", "--basis", dest="basis", default="def2tzvp", help="default: def2tzvp")
    parser.add_argument("--s", choices=["g", "o", "p"], default="g", help="g/o/p (default g)")
    add_job_args(parser)
    a = parser.parse_args()
    multi_spprep(a.method, a.basis, a.s, walltime=a.walltime, processors=a.cores, memory=a.mem, rerun=a.rerun)


if __name__ == "__main__":
    main_cli()
