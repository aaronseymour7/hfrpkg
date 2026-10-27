import glob
import os
import shutil

from hfrpkg.common import conformer_name, read_index, resolve_software, write_index
from hfrpkg.multi.grab_unique import numeric_key


def fill_logs(unique_folder="unique_files", subdir=""):
    """Copy each unique species' output into every *.mhfr folder that needs it.

    subdir="" handles optimisation logs (all conformers of a species); subdir="spec" handles single
    points (one per species, run on the lowest-energy conformer) and writes <n>.mhfr/spec/index.txt
    carrying the single-point software.
    """
    src_dir = os.path.join(unique_folder, subdir)
    index_path = os.path.join(src_dir, "index.txt")
    if not os.path.exists(index_path):
        print(f"Missing {index_path}")
        return

    unique = read_index(index_path)
    out_ext = resolve_software(unique.software)[2]
    by_inchi = {inchi: name for name, inchi, _ in unique.entries}

    for mhfr_dir in sorted(glob.glob("*.mhfr"), key=numeric_key):
        if not os.path.exists(os.path.join(mhfr_dir, "index.txt")):
            print(f"Skipping {mhfr_dir}, missing index.txt")
            continue
        index = read_index(mhfr_dir)
        dst_dir = os.path.join(mhfr_dir, subdir)
        os.makedirs(dst_dir, exist_ok=True)
        counts = dict(index.conformers)
        for name, inchi, _ in index.entries:
            if inchi not in by_inchi:
                print(f"[WARNING] InChI not found: {inchi}")
                continue
            uname = by_inchi[inchi]
            n = unique.conformers.get(uname, 1)
            pairs = [(uname, name)] if n <= 1 else [(conformer_name(uname, k), conformer_name(name, k)) for k in range(1, n + 1)]
            for src_stem, dst_stem in pairs:
                src = os.path.join(src_dir, src_stem + out_ext)
                if not os.path.exists(src):
                    print(f"[WARNING] Missing log file: {src}")
                    continue
                shutil.copyfile(src, os.path.join(dst_dir, dst_stem + out_ext))
            if not subdir:  # the unique set is authoritative for how many conformers were run
                if n > 1:
                    counts[name] = n
                else:
                    counts.pop(name, None)
        if subdir:
            write_index(dst_dir, index.level, unique.software, index.input_smiles, index.entries)
        else:
            write_index(mhfr_dir, index.level, index.software, index.input_smiles, index.entries,
                        conformers=counts, conformer_rule=unique.conformer_rule or index.conformer_rule)
