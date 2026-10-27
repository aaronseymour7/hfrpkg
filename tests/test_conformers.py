import os
import re

import numpy as np
import pytest
from rdkit import Chem

from hfrpkg import compute_folder as cf
from hfrpkg.common import conformer_name, read_index, write_index
from hfrpkg.conformers import ConformerSettings, generate_conformers, select_lowest
from hfrpkg.runner import run_reaction


def fake_energy(log):
    """c02 is the lowest conformer of every species; single-structure species get a fixed energy."""
    m = re.search(r"-c(\d+)\.", os.path.basename(log))
    e = -100.0 + 0.01 * abs(int(m.group(1)) - 2) if m else -100.0 - len(os.path.basename(log)) / 10
    return e, 0.05


def test_rule_matches_the_benchmark_scripts():
    s = ConformerSettings()
    assert [s.count(n) for n in (0, 1, 2, 3, 4)] == [5, 5, 10, 15, 20]  # max(5, 5 * n_rot)
    assert (s.seed, s.rmsd_cutoff, s.mmff_variant) == (42, 0.15, "MMFF94s")
    assert ConformerSettings(enabled=False).count(4) == 1


# species covering: no rotors, one/two/three rotors, heteroatoms, stereo, H2
INCHIS = [
    "InChI=1S/C4H10O/c1-2-3-4-5/h5H,2-4H2,1H3",          # 1-butanol
    "InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3",              # butane
    "InChI=1S/C5H12/c1-3-5-4-2/h3-5H2,1-2H3",            # pentane
    "InChI=1S/C4H10O/c1-3-5-4-2/h3-4H2,1-2H3",           # diethyl ether
    "InChI=1S/C4H8O/c1-3-4(2)5/h3H2,1-2H3",              # butan-2-one
    "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+",          # (E)-2-butene
    "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",                # benzene
    "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",                 # ethanol
    "InChI=1S/C3H6O2/c1-3(4)5-2/h1-2H3",                 # methyl acetate
    "InChI=1S/H2/h1H",                                   # dihydrogen
]


@pytest.mark.parametrize("inchi", INCHIS)
def test_generation_is_identical_to_the_original_script(inchi):
    from reference_makeconf import generate_conformers_from_inchi as reference

    ref_mol, ref_ids = reference(inchi)
    mol, ids = generate_conformers(inchi)
    assert len(ids) == len(ref_ids) and len(ids) >= 1
    for a, b in zip(ids, ref_ids):  # same conformers, same order, same coordinates
        assert np.allclose(mol.GetConformer(a).GetPositions(), ref_mol.GetConformer(b).GetPositions(), atol=1e-9)


def test_rotors_give_several_conformers_and_rigid_species_collapse():
    assert len(generate_conformers("InChI=1S/C5H12/c1-3-5-4-2/h3-5H2,1-2H3")[1]) > 1
    assert len(generate_conformers("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")[1]) == 1


def test_write_creates_every_conformer_and_records_counts(tmp_path):
    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCCCO", outfolder=str(out))
    idx = read_index(out)
    assert idx.conformer_rule and "EmbedMultipleConfs" in idx.conformer_rule and "seed=42" in idx.conformer_rule
    assert idx.conformers, "1-butanol and its fragments with rotors should have several conformers"
    for name, _, _ in idx.entries:
        n = idx.conformers.get(name, 1)
        stems = [name] if n == 1 else [conformer_name(name, k) for k in range(1, n + 1)]
        assert all((out / (s + ".com")).exists() for s in stems)
    # a re-run is reproducible: same files, same counts
    out2 = tmp_path / "2.mhfr"
    run_reaction("write", "isodesmic", "CCCCO", outfolder=str(out2))
    assert read_index(out2).conformers == idx.conformers
    assert (out / "R1_1-c01.com").read_text() == (out2 / "R1_1-c01.com").read_text() if idx.conformers.get("R1_1") else True


def test_no_conformers_flag_restores_single_structure_layout(tmp_path):
    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCCCO", outfolder=str(out), conformers=ConformerSettings(enabled=False))
    idx = read_index(out)
    assert not idx.conformers and idx.conformer_rule.startswith("single conformer")
    assert all((out / (n + ".com")).exists() for n, _, _ in idx.entries)


def _fake_run(folder, monkeypatch):
    idx = read_index(folder)
    for name, _, _ in idx.entries:
        n = idx.conformers.get(name, 1)
        for stem in ([name] if n == 1 else [conformer_name(name, k) for k in range(1, n + 1)]):
            (folder / (stem + ".log")).write_text("")
    monkeypatch.setattr(cf, "read_opt_energy", fake_energy)
    return idx


def test_compute_uses_lowest_energy_conformer_not_the_first(tmp_path, monkeypatch):
    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCCCO", outfolder=str(out))
    idx = _fake_run(out, monkeypatch)
    _, e, info = select_lowest(str(out), next(iter(idx.conformers)), next(iter(idx.conformers.values())), ".log")
    assert info["selected"].endswith("-c02.log") and e[0] == pytest.approx(-100.0)
    data = cf.compute_folder(str(out))
    assert data is not None and set(data["conformers"]) == set(idx.conformers)
    from hfrpkg.write_data import write_conformer_selection

    path = write_conformer_selection(data, str(out))
    text = open(path).read()
    assert "yes" in text and text.count("yes") == len(idx.conformers)


def test_missing_conformer_is_an_error_not_a_silent_partial_choice(tmp_path, monkeypatch):
    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCCCO", outfolder=str(out))
    idx = _fake_run(out, monkeypatch)
    name = next(iter(idx.conformers))
    os.remove(out / (conformer_name(name, 2) + ".log"))  # the lowest one never finished
    assert cf.compute_folder(str(out)) is None


def test_single_points_start_from_the_lowest_conformer(tmp_path, monkeypatch):
    import AaronTools.geometry as agm
    from hfrpkg import run_spec

    out = tmp_path / "1.mhfr"
    run_reaction("write", "isodesmic", "CCCCO", outfolder=str(out))
    idx = _fake_run(out, monkeypatch)
    used = {}

    class StubGeometry:
        def __init__(self, path):
            self.path = path

        def write(self, outfile, theory):
            used[os.path.basename(outfile)] = os.path.basename(self.path)
            open(outfile, "w").close()

    monkeypatch.setattr(agm, "Geometry", StubGeometry)
    run_spec.make_spec("m06-2x", "def2tzvp", "g", folder=str(out), submit=False)
    for name in idx.conformers:
        assert used[name + ".com"] == conformer_name(name, 2) + ".log"
    assert not read_index(out / "spec").conformers  # one single point per species


def test_multi_flow_carries_all_conformers(tmp_path, monkeypatch):
    from hfrpkg.multi.comp import main as comp
    from hfrpkg.multi.fill import fill_logs
    from hfrpkg.multi.grab_unique import grab_unique_coms
    from hfrpkg.multi.hfr_driver import main as prep

    monkeypatch.chdir(tmp_path)
    (tmp_path / "list.txt").write_text("2 CCCCO\n2 CCCC\n")
    prep("list.txt")
    grab_unique_coms()
    uidx = read_index("unique_files")
    assert uidx.conformers and uidx.conformer_rule
    for name, _, _ in uidx.entries:
        n = uidx.conformers.get(name, 1)
        for stem in ([name] if n == 1 else [conformer_name(name, k) for k in range(1, n + 1)]):
            assert os.path.exists(f"unique_files/{stem}.com")
            open(f"unique_files/{stem}.log", "w").close()
    monkeypatch.setattr(cf, "read_opt_energy", fake_energy)
    fill_logs("unique_files")
    for mhfr in ("1.mhfr", "2.mhfr"):
        idx = read_index(mhfr)
        for name, _, _ in idx.entries:
            n = idx.conformers.get(name, 1)
            for stem in ([name] if n == 1 else [conformer_name(name, k) for k in range(1, n + 1)]):
                assert os.path.exists(f"{mhfr}/{stem}.log")
    comp()
    assert os.path.getsize("enthalpies_summary.csv") > 100
