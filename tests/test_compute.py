import pytest

from hfrpkg import compute_folder as cf
from hfrpkg.common import (
    Index, get_Hf, inchi_from_smiles, parse_name, read_index, resolve_software, write_index,
)
from hfrpkg.read_optsum import read_optsum
from hfrpkg.write_data import write_single_reaction

ETOH = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"


def test_atct_data_is_packaged_and_found():
    assert get_Hf("InChI=1S/CH4/h1H4") == pytest.approx(-66.540)


def test_parse_name_and_software():
    assert parse_name("R3_2") == ("R", 3, 2)
    assert parse_name("junk") is None
    assert resolve_software("o") == ("ORCA", ".inp", ".out")
    assert resolve_software("Gaussian")[2] == ".log"
    with pytest.raises(ValueError):
        resolve_software("x")


def test_index_roundtrip(tmp_path):
    entries = [("R1_1", "InChI=1S/CH4/h1H4", "C"), ("P1_1", "InChI=1S/H2/h1H", "[H][H]")]
    write_index(tmp_path, "Isodesmic", "Psi4", "CCO", entries)
    idx = read_index(tmp_path)
    assert (idx.level, idx.software, idx.input_smiles) == ("Isodesmic", "Psi4", "CCO")
    assert idx.entries == entries
    assert idx.extensions == (".in", ".out")


def _index(order):
    ents = {
        "methane": ("InChI=1S/CH4/h1H4", "C"),
        "ethanol": (ETOH, "CCO"),
        "ethane": ("InChI=1S/C2H6/c1-2/h1-2H3", "CC"),
        "methanol": ("InChI=1S/CH4O/c1-2/h2H,1H3", "CO"),
    }
    names = {"methane": "R1_1", "ethanol": "R2_1", "ethane": "P1_1", "methanol": "P2_1"}
    return Index("Isodesmic", "Gaussian", "CCO",
                 [(names[k], *ents[k]) for k in order]), names


@pytest.mark.parametrize("order", [
    ["methane", "ethanol", "ethane", "methanol"],
    ["ethanol", "methane", "ethane", "methanol"],  # input not last
])
def test_input_is_identified_by_inchi_and_hf_is_consistent(order):
    index, names = _index(order)
    # R positions follow the order given, so renumber to match
    r = [e for e in index.entries if e[0].startswith("R")]
    index.entries = [(f"R{i}_1", inchi, s) for i, (_, inchi, s) in enumerate(r, 1)] + \
        [e for e in index.entries if e[0].startswith("P")]
    energies = {name: (-100.0 - i * 0.37, 0.05) for i, (name, _, _) in enumerate(index.entries)}
    data = cf.build_reaction(index, energies)

    assert data["input_inchi"] == ETOH
    assert data["reactants"][-1][2] == ETOH  # input reported last
    hf_p = sum(c * hf for c, _, _, hf, _, _ in data["products"])
    hf_r = sum(c * hf for c, _, _, hf, _, _ in data["reactants"][:-1])
    assert data["dft_hf"] == pytest.approx(hf_p - hf_r - data["reaction_H"], abs=1e-3)


def test_missing_energy_is_an_error_not_a_silent_partial_sum():
    index, _ = _index(["methane", "ethanol", "ethane", "methanol"])
    energies = {n: (-100.0, 0.05) for n, _, _ in index.entries}
    energies["P2_1"] = (None, None)
    with pytest.raises(cf.ComputeError, match="P2_1"):
        cf.build_reaction(index, energies)


def test_unknown_atct_gives_none_hf():
    index, _ = _index(["methane", "ethanol", "ethane", "methanol"])
    index.entries[0] = ("R1_1", "InChI=1S/NOT/IN/ATCT", "C")
    energies = {n: (-100.0, 0.05) for n, _, _ in index.entries}
    assert cf.build_reaction(index, energies)["dft_hf"] is None


def test_summary_roundtrip_with_missing_values(tmp_path):
    index, _ = _index(["methane", "ethanol", "ethane", "methanol"])
    index.entries[0] = ("R1_1", "InChI=1S/NOT/IN/ATCT", "C")
    energies = {n: (-100.0 - i, 0.05) for i, (n, _, _) in enumerate(index.entries)}
    data = cf.build_reaction(index, energies)
    row = write_single_reaction(data, tmp_path)
    assert row[3] is None  # no kcal value when dHf could not be determined
    back = read_optsum(tmp_path)
    assert back["dft_hf"] is None
    assert back["reactants"][0][3] is None
    assert back["reaction_H"] == pytest.approx(data["reaction_H"])


def test_inchi_from_smiles_matches_index_convention():
    assert inchi_from_smiles("CCO") == ETOH


def test_output_extensions_agree_with_aarontools():
    """Psi4 and ORCA write .out, Gaussian .log; a mismatch silently finds no output files."""
    from AaronTools.fileIO import expected_out_ext
    from hfrpkg.common import SOFTWARE

    for name, _, out_ext in SOFTWARE.values():
        assert out_ext == expected_out_ext(name)
