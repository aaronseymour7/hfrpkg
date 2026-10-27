"""Single points with a different program than the optimization (benchmark: Gaussian opt, Psi4 SP)."""
import os

import pytest

from hfrpkg import compute_folder as cf
from hfrpkg.common import read_index
from hfrpkg.conformers import ConformerSettings
from hfrpkg.runner import run_reaction


def test_gaussian_optimization_then_psi4_single_points(tmp_path, monkeypatch):
    import AaronTools.geometry as agm
    from hfrpkg import run_spec

    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCO", outfolder=str(out), conformers=ConformerSettings(enabled=False))
    idx = read_index(out)
    assert idx.software == "Gaussian"
    for name, _, _ in idx.entries:
        (out / (name + ".log")).write_text("")
    monkeypatch.setattr(cf, "read_opt_energy", lambda log: (-100.0, 0.05))

    class StubGeometry:
        def __init__(self, path):
            pass

        def write(self, outfile, theory):
            open(outfile, "w").close()

    monkeypatch.setattr(agm, "Geometry", StubGeometry)
    run_spec.make_spec("wb97x-v", "def2-tzvp", "p", folder=str(out), submit=False)

    spec = out / "spec"
    sp_idx = read_index(spec)
    assert sp_idx.software == "Psi4" and sp_idx.extensions == (".in", ".out")
    assert all((spec / (n + ".in")).exists() for n, _, _ in idx.entries)

    for name, _, _ in idx.entries:  # Psi4 output files (.out)
        (spec / (name + ".out")).write_text("")
    monkeypatch.setattr(cf, "read_sp_energy", lambda log, zpve: (-101.0 + zpve, zpve))
    data = cf.compute_spec_folder(str(out))
    assert data is not None and data["input_inchi"].startswith("InChI=1S/C2H6O/")
