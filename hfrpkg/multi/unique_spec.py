import glob
import os
import sys

from hfrpkg.common import read_index, resolve_software, write_index
from hfrpkg.conformers import select_lowest
from hfrpkg.multi.run_unique import run_jobs


def make_spec(method, basis, extension, folder="unique_files", submit=True, **job_kwargs):
    """Single-point inputs for the unique molecules, built from their optimised logs."""
    from AaronTools.geometry import Geometry
    from AaronTools.theory import Theory
    from AaronTools.theory.job_types import SinglePointJob

    unique = read_index(folder)
    _, opt_out = unique.extensions
    sp_software, sp_in, _ = resolve_software(extension)

    if not glob.glob(os.path.join(folder, "*" + opt_out)):
        print(f"No {opt_out} files found in {folder}.")
        sys.exit(1)

    spec_dir = os.path.join(folder, "spec")
    os.makedirs(spec_dir, exist_ok=True)
    write_index(spec_dir, unique.level, sp_software, None, unique.entries)

    theory = Theory(method=method, basis=basis, job_type=SinglePointJob())
    for name, _, _ in unique.entries:
        try:
            log, _, _ = select_lowest(folder, name, unique.conformers.get(name, 1), opt_out)
            Geometry(log).write(outfile=os.path.join(spec_dir, name + sp_in), theory=theory)
        except Exception as exc:
            print(f"Error processing {name}: {exc}")

    if submit:
        run_jobs(spec_dir, **job_kwargs)
