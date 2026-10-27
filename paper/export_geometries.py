"""Cartesian coordinates and absolute energies of the selected conformer of every species, for the SI.

    python export_geometries.py RUN_DIR [RUN_DIR ...] [--out si_structures]

Each RUN_DIR is a multi-run directory (N.mhfr folders after multi-compute, and N.mhfr/spec after
multi-spcomp).  For every species, the lowest-energy conformer (E + ZPVE) is taken, as in the
thermochemistry.  Writes:

  <out>.xyz  one block per species and run: atom count, a comment line with the run directory, name,
             SMILES, InChI, energies (hartree) and number of imaginary frequencies, then the coordinates
  <out>.csv  run_dir, smiles, inchi, conformer, E_hartree, ZPVE_hartree, E_plus_ZPVE_hartree,
             n_imaginary, E_sp_hartree (single point on that geometry, blank if none)

Species shared by several reactions are written once per run directory, unless reactions computed
the same species separately and got different energies; then each distinct structure is written
(the conformer field, N.mhfr/file, says which reaction folder it came from).  Any species with an
imaginary frequency is also printed, so that the "no imaginary frequencies" statement can be checked.
"""
import argparse
import csv
import glob
import os

from hfrpkg.common import read_index, resolve_software
from hfrpkg.conformers import select_lowest


def reader(path):
    from AaronTools.fileIO import FileReader

    return FileReader(path, just_geom=False)


def species_rows(run_dir):
    seen = {}
    for folder in sorted(glob.glob(os.path.join(run_dir, "*.mhfr")), key=lambda p: int(os.path.basename(p).split(".")[0])):
        index = read_index(folder)
        _, out_ext = index.extensions
        spec_dir = os.path.join(folder, "spec")
        sp_ext = None
        if os.path.exists(os.path.join(spec_dir, "index.txt")):
            sp_ext = resolve_software(read_index(spec_dir).software)[2]
        for name, inchi, smiles in index.entries:
            try:
                log, (e_zpve, zpve), _ = select_lowest(folder, name, index.conformers.get(name, 1), out_ext)
            except ValueError as exc:
                print(f"[WARN] {folder}: {exc}")
                continue
            key = (inchi, round(e_zpve, 6))
            if key in seen:
                continue
            r = reader(log)
            imag = r["frequency"].imaginary_frequencies if "frequency" in r.keys() else None
            n_imag = None if imag is None else len(imag)
            e_sp = None
            if sp_ext and os.path.exists(os.path.join(spec_dir, name + sp_ext)):
                sp = reader(os.path.join(spec_dir, name + sp_ext))
                e_sp = sp["energy"] if "energy" in sp.keys() else None
            seen[key] = dict(run_dir=run_dir, smiles=smiles, inchi=inchi,
                             conformer=os.path.basename(folder) + "/" + os.path.basename(log),
                             E=e_zpve - zpve, ZPVE=zpve, E_ZPVE=e_zpve, n_imag=n_imag, imag=imag or [], E_sp=e_sp,
                             atoms=[(a.element, *a.coords) for a in r["atoms"]])
    return list(seen.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--out", default="si_structures")
    args = ap.parse_args()
    rows = [row for d in args.run_dirs for row in species_rows(d)]
    with open(args.out + ".xyz", "w") as fh:
        for r in rows:
            fh.write(f"{len(r['atoms'])}\n")
            sp = f" E_sp={r['E_sp']:.8f}" if r["E_sp"] is not None else ""
            fh.write(f"{r['run_dir']} {r['conformer']} {r['smiles']} {r['inchi']} E={r['E']:.8f} ZPVE={r['ZPVE']:.6f}"
                     f" E+ZPVE={r['E_ZPVE']:.8f}{sp} n_imag={r['n_imag']}\n")
            for el, x, y, z in r["atoms"]:
                fh.write(f"{el:<2} {x:14.8f} {y:14.8f} {z:14.8f}\n")
    with open(args.out + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_dir", "smiles", "inchi", "conformer", "E_hartree", "ZPVE_hartree", "E_plus_ZPVE_hartree", "n_imaginary", "E_sp_hartree"])
        for r in rows:
            w.writerow([r["run_dir"], r["smiles"], r["inchi"], r["conformer"], f"{r['E']:.8f}", f"{r['ZPVE']:.6f}",
                        f"{r['E_ZPVE']:.8f}", r["n_imag"], "" if r["E_sp"] is None else f"{r['E_sp']:.8f}"])
    multi = {}
    for r in rows:
        multi.setdefault((r["run_dir"], r["inchi"]), []).append(r)
    multi = {k: v for k, v in multi.items() if len(v) > 1}
    print(f"species with different energies in different reactions: {len(multi)}")
    for (d, _), v in multi.items():
        spread = (max(x["E_ZPVE"] for x in v) - min(x["E_ZPVE"] for x in v)) * 627.5095
        print(f"  {d} {v[0]['smiles']}: " + ", ".join(x["conformer"] for x in v) + f" (spread {spread:.2f} kcal/mol)")
    imag = [r for r in rows if r["n_imag"]]
    print(f"{len(rows)} structures written to {args.out}.xyz / .csv; with imaginary frequencies: {len(imag)}")
    for r in imag:
        print(f"  {r['run_dir']} {r['smiles']} {r['conformer']}: {r['n_imag']} imaginary ("
              + ", ".join(f"{abs(v):.1f}i" for v in r["imag"]) + " cm-1)")


if __name__ == "__main__":
    main()
