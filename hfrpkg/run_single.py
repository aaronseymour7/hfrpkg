import argparse
import glob
import os
import sys

from hfrpkg.common import read_index, working_dir


def run_jobs(folder=".", walltime=12, processors=8, memory=12, rerun=False):
    """Submit every input file in `folder`; skip ones whose output already exists."""
    from AaronTools.job_control import SubmitProcess

    try:
        in_ext, out_ext = read_index(folder).extensions
    except (OSError, ValueError) as exc:
        print(f"Could not determine software from {folder}/index.txt: {exc}")
        sys.exit(1)

    with working_dir(folder):
        inputs = sorted(glob.glob("*" + in_ext))
        if not inputs:
            print(f"No {in_ext} files found in {folder}.")
            sys.exit(1)
        for f in inputs:
            if not rerun and os.path.exists(os.path.splitext(f)[0] + out_ext):
                print(f"skip {f}: output exists (use --rerun to force)")
                continue
            try:
                SubmitProcess(f, walltime, processors, memory).submit()
            except Exception as exc:
                print(f"failed to submit {f}: {exc}")


def add_job_args(parser):
    parser.add_argument("--walltime", type=int, default=12, help="hours (default 12)")
    parser.add_argument("--cores", type=int, default=8, help="processors (default 8)")
    parser.add_argument("--mem", type=int, default=12, help="GB (default 12)")
    parser.add_argument("--rerun", action="store_true", help="resubmit even if output exists")


def main_cli():
    parser = argparse.ArgumentParser(description="Submit the jobs in a reaction folder")
    parser.add_argument("folder", nargs="?", default=".")
    add_job_args(parser)
    a = parser.parse_args()
    run_jobs(a.folder, a.walltime, a.cores, a.mem, a.rerun)


if __name__ == "__main__":
    main_cli()
