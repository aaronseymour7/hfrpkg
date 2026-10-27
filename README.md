# hfrpkg [![DOI](https://zenodo.org/badge/1435294364.svg)](https://doi.org/10.5281/zenodo.23683594)


Generate balanced **isogyric, isodesmic, hypohomodesmotic and homodesmotic** reactions for a molecule, run them
with Gaussian / ORCA / Psi4, and turn the results into an enthalpy of formation referenced to
[ATcT](https://atct.anl.gov) values. Built on RDKit, PuLP and AaronTools.

## Install

```bash
pip install -e ".[test]"      # RDKit, PuLP (<4, bundles CBC), AaronTools, pytest
pytest                        # no QM software needed
```

## One molecule

```bash
hfr view  isodesmic CCO                       # print the balanced reaction
hfr count homodesmotic c1ccccc1               # show the bond/atom counts being balanced
hfr write homodesmotic CCCCO --outfolder 1.mhfr --m B3LYP --b 6-31G --s g

folder-run 1.mhfr                             # submit (skips jobs whose output exists; --rerun to force)
folder-compute 1.mhfr                         # opt_summary.txt + enthalpy_summary.csv

folder-sprun 1.mhfr --m wb97x-v --b def2-tzvp --s p   # Psi4 single points on the optimised geometries -> 1.mhfr/spec
folder-spcomp --f 1.mhfr                      # spec/sp_summary.txt
```

Job resources: `--walltime` (h), `--cores`, `--mem` (GB) on every submit command.

## Many molecules

`list.txt`: one `<level> <SMILES>` per line (1 isogyric, 2 isodesmic, 3 hypohomodesmotic, 4 homodesmotic).

```bash
multi-prep list.txt          # writes N.mhfr/ for each line (N = line number)
multi-run                    # one job per *unique* molecule -> unique_files/
multi-compute                # copy logs back into each N.mhfr, write Reactions_summary.txt, enthalpies_summary.csv
multi-sprun                  # single points for the unique molecules -> unique_files/spec
multi-spcomp                 # sp_reaction_summaries.txt, sp_enthalpies_summary.csv
```

## Coming from HFR_PACKAGE

The commands and their arguments are those of the original
[HFR_PACKAGE](https://github.com/aaronseymour7/HFR_PACKAGE); every original command line runs unchanged here.
The two conformer scripts that were run by hand between the original steps are built in:

| Original workflow | Here |
|---|---|
| `multi-prep list.txt --m M --b B` | same (also writes the conformers) |
| `python makeConf.py` (in the run directory: conformers of `unique_files/`, submitted 24 h / 12 cores / 24 GB) | `multi-run --walltime 24 --cores 12 --mem 24` |
| `python extractConfs.py` (lowest E + ZPVE copied into `unique_files/`) | done by `multi-compute` |
| `multi-compute`, `multi-sprun --m M --b B --s p`, `multi-spcomp` | same |

Psi4 outputs are `.out` (the original expected `.dat`).

## Conformers

Each species is optimized from several starting conformers; the **lowest-energy conformer** (E + ZPVE from
its frequency job) is used for the thermochemistry and as the geometry for single points. The procedure is
a port of the `makeConf.py` / `extractConfs.py` scripts used for the published benchmark, and
`tests/test_conformers.py` checks it conformer-by-conformer against a verbatim copy of the original:

1. Start from the species' InChI, add hydrogens.
2. Embed `n = max(5, 5 × n_rotatable_bonds)` conformers with RDKit `EmbedMultipleConfs` (seed 42).
3. Prune by best-fit all-atom RMSD (0.15 Å), keeping the first of each cluster.
4. Minimize with MMFF94s and order by force-field energy.
5. Prune again (0.15 Å) in that order, so the lowest-energy member of each cluster survives.

Species left with one conformer keep `R1_2.com`; others get `R1_2-c01.com … R1_2-cNN.com` (c01 = lowest
force-field energy). `index.txt` records the count and the settings, including the RDKit version (the
embedded conformers depend on it). Options: `--no-conformers`, `--min-conformers`, `--per-rotor`,
`--conf-rmsd`, `--conf-seed`.

**RDKit version.** The embedded conformers depend on the RDKit version. The benchmark used
**RDKit 2024.03.3** (`pip install rdkit==2024.3.3`, Python ≤ 3.12). A newer RDKit (2026.03.6) gives a
different conformer set for 13 of the 72 benchmark molecules (182 vs 173 conformer jobs) and a different
lowest force-field energy for 3. `paper/data/conformer_counts_rdkit2024.03.3.csv` lists the expected
conformer count and MMFF94s energies of each benchmark molecule, and CI checks them under that version.

- `folder-run` / `multi-run` submit every conformer; `folder-compute` / `multi-compute` pick the lowest and
  write `conformer_selection.csv` (all conformer energies and the one used). An unfinished conformer is an
  error, because choosing among a partial set could silently pick the wrong one.
- `folder-sprun` / `multi-sprun` start from the selected conformer, one single point per species.

## Software used for the published benchmark

| Component | Version |
|---|---|
| Geometry optimization and frequencies | Gaussian 16, Revision C.01 |
| Single-point energies | Psi4 1.9.1 |
| RDKit (conformers) | 2024.03.3 |
| PuLP / CBC (reaction ILP) | 2.8.0 / 2.10.3 |
| ATcT reference values | `hfrpkg/data/ATcT_lib.txt` (ATcT version 1.156; uncertainties are 95% confidence intervals) |

The reactions generated for the benchmark were re-checked under PuLP 2.8.0 and RDKit 2024.03.3
(`paper/objective_check.py`): the isodesmic and homodesmotic reactions are independent of the solver
formulation; see the manuscript, Table S4, for the four hypohomodesmotic exceptions.

## Reproducing the paper

Everything in `paper/` regenerates from the data in the repository (run from `paper/`, with the package installed):

```bash
python make_figures.py          # figures/*.png (600 dpi) and *.pdf, from data/benchmark_errors.csv and data/scaling.csv
python paper_statistics.py      # every statistic quoted in the text
python scaling_demo.py          # data/scaling.csv (polyenes, PAHs, C60)
python objective_check.py PATH/TO/legacy/hfrpkg/core.py   # data/objective_check.csv
python atct_uncertainty.py      # data/atct_uncertainty.csv: ATcT uncertainties of the targets and propagated reference uncertainties
python wildcard_libraries.py    # data/wildcard_libraries.csv: the fragment libraries of each level, extracted from core.py
python feasibility_scan.py      # data/feasibility.csv: applicability survey on 229 ATcT species (identity-only vs no reaction)
python selection_check.py       # data/selection_check.csv: which survey species meet the benchmark's selection criteria
python manuscript_inputs.py     # data/stats.json, mols.json, hy_ambiguous.json for the manuscript builder
python export_geometries.py RUN_DIR ...   # SI coordinates and absolute energies of the selected conformers (needs the output files)
python reopt_saddles.py RUN_DIR ...     # restart optimizations that ended on a saddle point, displaced along the imaginary modes
node build_manuscript.js        # Seymour_Wheeler_Homodesmotic_manuscript.docx  (needs the npm package `docx`)
```

`export_data.py` rebuilds `data/benchmark_errors.csv` from the original result workbooks (not in the repository).

## Method notes

- Energies are **E + ZPVE** (0 K); the ATcT column used is the 0 K ΔfH°, so the two are consistent.
- Single points use the single-point electronic energy plus the ZPVE from the optimisation.
- The input molecule is found by InChI, not by file order. Missing energies abort that folder
  (a partial sum would silently give a wrong reaction enthalpy). A species absent from ATcT leaves
  ΔHf blank rather than guessing.
- The solver minimises the number of species and never returns the identity reaction
  (input → input). For small molecules the homodesmotic classes can therefore be *infeasible* — that is
  reported as such.

## index.txt format

```
Level:<TAB>Isodesmic<TAB>Software:<TAB>Gaussian
Input SMILES:<TAB>CCO
Filename<TAB>InChI<TAB>SMILES
R1_1 ...        # R = reactant, P = product, _N = stoichiometric coefficient
```
