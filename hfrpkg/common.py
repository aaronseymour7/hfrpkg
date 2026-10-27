"""Shared helpers: software tables, index.txt I/O, ATcT lookup, solver choice.

Everything that used to be copy-pasted across compute_folder / spec_compute /
fill / unique_spec lives here so the pieces cannot drift apart.
"""

import functools
import importlib.resources as resources
import os
import re
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

HARTREE_TO_KJ = 2625.4996
KJ_PER_KCAL = 4.184

# code -> (name, input extension, output extension)
SOFTWARE = {
    "g": ("Gaussian", ".com", ".log"),
    "o": ("ORCA", ".inp", ".out"),
    "p": ("Psi4", ".in", ".out"),
}
_BY_NAME = {v[0].lower(): v for v in SOFTWARE.values()}

NAME_RE = re.compile(r"^([PR])(\d+)_(\d+)$")


def resolve_software(key):
    """Accept 'g'/'o'/'p' or 'Gaussian'/'ORCA'/'Psi4'; return (name, in_ext, out_ext)."""
    if key is None:
        key = "g"
    k = str(key).strip().lower()
    if k in SOFTWARE:
        return SOFTWARE[k]
    if k in _BY_NAME:
        return _BY_NAME[k]
    raise ValueError(f"Unknown software {key!r}; use one of g, o, p")


def parse_name(name):
    """'R2_3' -> ('R', 2, 3): kind, position, stoichiometric coefficient."""
    m = NAME_RE.match(name)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


@contextmanager
def working_dir(path):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


# ---------------------------------------------------------------- index.txt
@dataclass
class Index:
    level: Optional[str] = None
    software: Optional[str] = None
    input_smiles: Optional[str] = None
    entries: List[Tuple[str, str, str]] = field(default_factory=list)  # (name, inchi, smiles)
    conformers: dict = field(default_factory=dict)  # name -> number of conformers (absent/0: one structure)
    conformer_rule: Optional[str] = None  # how the conformer set was generated, for reproducibility

    @property
    def extensions(self):
        """(input extension, output extension) for this index's software."""
        _, i, o = resolve_software(self.software)
        return i, o

    def inchi_of(self, name):
        for n, inchi, _ in self.entries:
            if n == name:
                return inchi
        return None

    def smiles_of(self, name):
        for n, _, smi in self.entries:
            if n == name:
                return smi
        return None


def read_index(path):
    """Parse an index.txt (file or the folder holding it)."""
    if os.path.isdir(path):
        path = os.path.join(path, "index.txt")
    idx = Index()
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            head = parts[0].strip()
            if head == "Level:":
                idx.level = parts[1] if len(parts) > 1 else None
                if "Software:" in parts:
                    i = parts.index("Software:") + 1
                    if i < len(parts):
                        idx.software = parts[i].strip()
            elif head.startswith("Software:"):  # legacy "Software: X" header
                idx.software = (parts[1] if len(parts) > 1 else head.split(":", 1)[1]).strip()
            elif head == "Input SMILES:":
                idx.input_smiles = parts[1] if len(parts) > 1 else None
            elif head == "Conformers:":
                idx.conformer_rule = parts[1] if len(parts) > 1 else None
            elif head == "Filename":
                continue
            elif len(parts) >= 2:
                idx.entries.append((parts[0], parts[1], parts[2] if len(parts) > 2 else ""))
                if len(parts) > 3 and parts[3].strip().isdigit() and int(parts[3]) > 1:
                    idx.conformers[parts[0]] = int(parts[3])
    return idx


def conformer_name(name, k):
    """File stem of conformer k (1-based) of species `name`, e.g. R2_1-c03."""
    return f"{name}-c{k:02d}"


def write_index(path, level, software, input_smiles, entries, conformers=None, conformer_rule=None):
    """entries: (name, inchi, smiles); conformers: {name: n} for species with more than one structure."""
    if os.path.isdir(path):
        path = os.path.join(path, "index.txt")
    conformers = conformers or {}
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"Level:\t{level}\tSoftware:\t{software}\n")
        if input_smiles is not None:
            fh.write(f"Input SMILES:\t{input_smiles}\n")
        if conformer_rule:
            fh.write(f"Conformers:\t{conformer_rule}\n")
        fh.write("Filename\tInChI\tSMILES\tConformers\n")
        for name, inchi, smiles in entries:
            fh.write(f"{name}\t{inchi}\t{smiles}\t{conformers.get(name, 1)}\n")


def inchi_from_smiles(smiles):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchi(Chem.AddHs(mol))


# --------------------------------------------------------------------- ATcT
@functools.lru_cache(maxsize=1)
def _atct_table():
    """InChI -> dHf(0 K, kJ/mol). Column 6 of ATcT_lib.txt is the 0 K value, which
    is the right partner for E + ZPVE energies."""
    text = resources.files("hfrpkg.data").joinpath("ATcT_lib.txt").read_text(encoding="utf-8")
    table = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 6:
            try:
                table[parts[3]] = float(parts[5])
            except ValueError:
                continue
    return table


def get_Hf(inchi):
    return _atct_table().get(inchi)


# ------------------------------------------------------------------- solver
def make_solver():
    """Prefer the bundled CBC; fall back to HiGHS if only that is present (pip install "pulp<4")."""
    import pulp

    if hasattr(pulp, "PULP_CBC_CMD"):
        cbc = pulp.PULP_CBC_CMD(msg=False)
        if cbc.available():
            return cbc
    for name in ("HiGHS", "HiGHS_CMD", "COIN_CMD"):
        if hasattr(pulp, name):
            solver = getattr(pulp, name)(msg=False)
            if solver.available():
                return solver
    raise RuntimeError("No MILP solver available. Install one with: pip install 'pulp<4'")
