"""Restart every optimization in a multi-run that ended on a saddle point (imaginary frequencies).

    python reopt_saddles.py RUN_DIR [RUN_DIR ...] [--step 0.15]

For each Gaussian output in RUN_DIR/unique_files with imaginary frequencies, the geometry is displaced
along the imaginary normal modes (the largest atomic displacement of each mode is --step Å, default 0.15)
and a new input is written in place of the old one, with "Opt" changed to "Opt=(CalcFC,Tight)" and
"NoSymm" added so that the optimizer cannot stay on the symmetric saddle point. Method, basis set and
integration grid are unchanged. The copies of the input in the N.mhfr folders are replaced too, so that
multi-run (which recopies them into unique_files) submits the new one. The old output goes to
unique_files/saddle_backup/opt/ and the single point of that species to unique_files/saddle_backup/spec/.
Afterwards run, in RUN_DIR:

    multi-run && (when finished) multi-compute && multi-sprun ... && (when finished) multi-spcomp
"""
import argparse
import glob
import os
import re
import shutil

import numpy as np

COORD = re.compile(r"^\s*([A-Z][a-z]?)\s+(-?\d+\.\d*)\s+(-?\d+\.\d*)\s+(-?\d+\.\d*)\s*$")


def displaced(log, step):
    from AaronTools.fileIO import FileReader

    r = FileReader(log, just_geom=False)
    freq = r["frequency"] if "frequency" in r.keys() else None
    if freq is None or not freq.imaginary_frequencies:
        return None, []
    xyz = np.array([a.coords for a in r["atoms"]], dtype=float)
    for mode in freq.data:
        if mode.frequency < 0:
            v = np.asarray(mode.vector, dtype=float)
            xyz = xyz + step * v / np.abs(np.linalg.norm(v, axis=1)).max()
    return xyz, freq.imaginary_frequencies


def rewrite_input(com, xyz):
    lines = open(com).read().splitlines()
    out, k = [], 0
    for line in lines:
        if line.startswith("#"):
            line = re.sub(r"\bopt\b(=\S+)?", "Opt=(CalcFC,Tight)", line, flags=re.I)
            if not re.search(r"\bnosymm\b", line, re.I):
                line = line.rstrip() + " NoSymm"
        else:
            m = COORD.match(line)
            if m and k < len(xyz):
                line = f"{m.group(1):<2} {xyz[k][0]:14.8f} {xyz[k][1]:14.8f} {xyz[k][2]:14.8f}"
                k += 1
        out.append(line)
    if k != len(xyz):
        raise ValueError(f"{com}: {k} coordinate lines but {len(xyz)} atoms")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--step", type=float, default=0.15)
    args = ap.parse_args()
    for run_dir in args.run_dirs:
        unique = os.path.join(run_dir, "unique_files")
        backup = os.path.join(unique, "saddle_backup")
        for log in sorted(glob.glob(os.path.join(unique, "*.log"))):
            com = log[:-4] + ".com"
            if not os.path.exists(com):
                print(f"[WARN] no input for {log}")
                continue
            xyz, imag = displaced(log, args.step)
            if xyz is None:
                continue
            old = open(com).read()
            text = rewrite_input(com, xyz)
            for sub in ("opt", "spec"):
                os.makedirs(os.path.join(backup, sub), exist_ok=True)
            name = os.path.basename(log)
            shutil.move(log, os.path.join(backup, "opt", name))
            shutil.copy(com, os.path.join(backup, "opt", name[:-4] + ".com"))
            copies = [c for c in glob.glob(os.path.join(run_dir, "*.mhfr", "*.com")) if open(c).read() == old]
            for c in [com] + copies:
                open(c, "w").write(text)
            species = re.sub(r"-c\d+$", "", name[:-4])
            for sp in glob.glob(os.path.join(unique, "spec", species + ".*")):
                shutil.move(sp, os.path.join(backup, "spec", os.path.basename(sp)))
            print(f"{log}: {len(imag)} imaginary (" + ", ".join(f"{abs(f):.1f}i" for f in imag)
                  + f" cm-1) -> new input written ({len(copies)} N.mhfr copies updated)")


if __name__ == "__main__":
    main()
