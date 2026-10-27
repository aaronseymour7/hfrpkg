import os

import pytest

from hfrpkg import compute_folder as cf
from hfrpkg.common import conformer_name, read_index
from hfrpkg.conformers import ConformerSettings
from hfrpkg.runner import run_reaction


@pytest.mark.parametrize("software,ext", [("g", ".com"), ("o", ".inp"), ("p", ".in")])
def test_write_creates_inputs_and_index(tmp_path, software, ext):
    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCO", outfolder=str(out), extension=software)
    idx = read_index(out)
    assert idx.input_smiles == "CCO" and idx.level == "Isodesmic"
    assert idx.entries
    for name, _, _ in idx.entries:
        n = idx.conformers.get(name, 1)  # species with several conformers are written as name-cNN
        stems = [name] if n == 1 else [conformer_name(name, k) for k in range(1, n + 1)]
        assert all((out / (s + ext)).exists() for s in stems)
    # input is the last reactant
    last_r = [e for e in idx.entries if e[0].startswith("R")][-1]
    assert last_r[2] == "CCO" or last_r[1].startswith("InChI=1S/C2H6O/")


def test_full_folder_roundtrip_with_fake_logs(tmp_path, monkeypatch):
    out = tmp_path / "1.mhfr"
    run_reaction("write", "isogyric", "CCO", outfolder=str(out), conformers=ConformerSettings(enabled=False))  # H2 listed after input in the solver
    idx = read_index(out)
    for name, _, _ in idx.entries:
        (out / (name + ".log")).write_text("")
    monkeypatch.setattr(cf, "read_opt_energy", lambda log: (-100.0 - len(os.path.basename(log)) / 10, 0.04))
    data = cf.compute_folder(str(out))
    assert data is not None
    assert data["input_smiles"] == "CCO"
    assert data["reactants"][-1][2].startswith("InChI=1S/C2H6O/")
    assert data["dft_hf"] is not None
