#!/usr/bin/env python3
import argparse

from hfrpkg.multi.grab_unique import grab_unique_coms
from hfrpkg.multi.run_unique import run_jobs
from hfrpkg.run_single import add_job_args


def main(**job_kwargs):
    grab_unique_coms()
    run_jobs(**job_kwargs)


def main_cli():
    parser = argparse.ArgumentParser(description="Submit one job per unique molecule across all *.mhfr folders")
    add_job_args(parser)
    a = parser.parse_args()
    main(walltime=a.walltime, processors=a.cores, memory=a.mem, rerun=a.rerun)


if __name__ == "__main__":
    main_cli()
