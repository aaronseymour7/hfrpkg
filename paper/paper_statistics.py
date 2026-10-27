"""Every statistic quoted in the manuscript text, computed from data/benchmark_errors.csv.

Run:  python paper_statistics.py
"""
import csv
import os
import math
from math import comb

import numpy as np

rows = list(csv.DictReader(open("data/benchmark_errors.csv", encoding="utf-8")))
mols = []
for r in rows:
    if r["inchi"] not in mols:
        mols.append(r["inchi"])
LV = sorted({r["geometry_level"] for r in rows})
cls = {r["inchi"]: r["class"] for r in rows}
S = ["ID", "HY", "HO"]
E = {(st, lv, s): np.array([float(next(r[s] for r in rows if r["inchi"] == m and r["stage"] == st and r["geometry_level"] == lv)) for m in mols])
     for st in ("opt", "sp") for lv in LV for s in S}
strained = np.array([cls[m] == "strained" for m in mols])
SUB = {"unstrained": ~strained, "strained": strained, "all": np.ones(len(mols), bool)}


def pooled(st, s, sub):
    return np.concatenate([E[(st, lv, s)][SUB[sub]] for lv in LV])


rmse = lambda v: math.sqrt(np.mean(v ** 2))
rng = np.random.default_rng(7)


def boot_diff(st, a, b, sub, n=4000):
    idx0 = np.where(SUB[sub])[0]
    out = []
    for _ in range(n):
        idx = rng.choice(idx0, len(idx0))
        va = np.concatenate([E[(st, lv, a)][idx] for lv in LV]); vb = np.concatenate([E[(st, lv, b)][idx] for lv in LV])
        out.append(rmse(va) - rmse(vb))
    return np.percentile(out, [2.5, 97.5])


def sign_p(k, n):
    k = max(k, n - k)
    return min(1.0, 2 * sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n)


print("== Tables 1/2: pooled RMSE / MAE / ME / STD(n-1)")
for st in ("opt", "sp"):
    for sub in ("unstrained", "strained", "all"):
        for s in S:
            v = pooled(st, s, sub)
            print(f"{st:3s} {sub:10s} {s}  N={len(v):3d} RMSE={rmse(v):.3f} MAE={np.mean(np.abs(v)):.3f} ME={v.mean():+.3f} STD={v.std(ddof=1):.3f}")
print("\n== within 1 kcal/mol (all, SPR) and max |error|")
for s in S:
    v = pooled("sp", s, "all")
    print(s, f"{100 * np.mean(np.abs(v) <= 1):.1f}%  max {np.abs(v).max():.2f}")
print("\n== HY vs HO")
for st in ("opt", "sp"):
    same = sum(abs(E[(st, lv, "HY")][i] - E[(st, lv, "HO")][i]) < 1e-9 for lv in LV for i in range(len(mols)))
    diffmol = {i for lv in LV for i in range(len(mols)) if abs(E[(st, lv, "HY")][i] - E[(st, lv, "HO")][i]) > 1e-9}
    print(f"{st}: identical {same}/{sum(r['stage'] == st for r in rows)}; molecules where they differ: {len(diffmol)} (acyclic {sum(not strained[i] for i in diffmol)}, strained {sum(strained[i] for i in diffmol)})")
    for sub in ("unstrained", "strained", "all"):
        pairs = [(abs(E[(st, lv, "HY")][i]), abs(E[(st, lv, "HO")][i])) for lv in LV for i in np.where(SUB[sub])[0] if abs(E[(st, lv, "HY")][i] - E[(st, lv, "HO")][i]) > 1e-9]
        b = sum(o < h for h, o in pairs)
        lo, hi = boot_diff(st, "HY", "HO", sub)
        print(f"  {sub:10s} differing {len(pairs)}, HO better {b}, sign-test p={sign_p(b, len(pairs)):.3f}; RMSE(HY)-RMSE(HO)={rmse(pooled(st,'HY',sub))-rmse(pooled(st,'HO',sub)):.3f} [{lo:.3f}, {hi:.3f}]")
print("\n== ID - HO (SPR)")
for sub in ("unstrained", "strained", "all"):
    lo, hi = boot_diff("sp", "ID", "HO", sub)
    print(f"  {sub:10s} {rmse(pooled('sp','ID',sub))-rmse(pooled('sp','HO',sub)):.3f} [{lo:.3f}, {hi:.3f}]")
print("\n== geometry spread of HO error (max-min over 3 levels)")
for st in ("opt", "sp"):
    st_ = np.array([E[(st, lv, "HO")] for lv in LV]); sp = st_.max(0) - st_.min(0)
    print(f"  {st}: mean {sp.mean():.2f}  max {sp.max():.1f}")

if os.path.exists("data/atct_uncertainty.csv"):
    u = list(csv.DictReader(open("data/atct_uncertainty.csv", encoding="utf-8")))
    col = lambda c: np.array([float(r[c]) for r in u])
    print("\n== ATcT uncertainty (kcal/mol; target as listed, reference species propagated in quadrature)")
    for c in ("sigma_target_kcal", "sigma_refs_ID", "sigma_refs_HY", "sigma_refs_HO"):
        v = col(c)
        print(f"  {c:18s} median {np.median(v):.2f}  mean {v.mean():.2f}  max {v.max():.2f}  RMS {math.sqrt(np.mean(v ** 2)):.2f}")
    comb = math.sqrt(np.mean(col("sigma_target_kcal") ** 2 + col("sigma_refs_HO") ** 2))
    print(f"  combined RMS floor, target + HO references: {comb:.2f}")
    print("  molecules with listed uncertainty > 0.5:", [r["smiles"] for r in u if float(r["sigma_target_kcal"]) > 0.5])

# ---- molecule-level inference (the unit of independent observation is the molecule, not the reaction)
from math import comb as _comb


def _sign_p(k, n):
    k = max(k, n - k)
    return min(1.0, 2 * sum(_comb(n, i) for i in range(k, n + 1)) / 2 ** n)


def _mean_abs(stage, m, s):
    return float(np.mean([abs(E[(stage, lv, s)][m]) for lv in LV]))


print("\n== molecule-level sign tests (mean |error| over the three levels; molecules where the schemes differ)")
for a, b in (("HY", "HO"), ("ID", "HO")):
    for st_ in ("opt", "sp"):
        if a == "ID" and st_ == "opt":
            continue
        for sub in ("unstrained", "strained", "all"):
            idx = np.where(SUB[sub])[0]
            d = [_mean_abs(st_, i, a) - _mean_abs(st_, i, b) for i in idx]
            diff = [x for x in d if abs(x) > 1e-9]
            k = sum(x > 0 for x in diff)
            print(f"  {a} vs {b}, {st_:3s} {sub:10s}: {b} better for {k:2d} of {len(diff):2d} molecules that differ, p={_sign_p(k, len(diff)):.3f}")
same = sum(all(abs(E[("sp", lv, "HY")][i] - E[("sp", lv, "HO")][i]) < 1e-9 for lv in LV) for i in range(len(mols)))
print(f"  molecules with identical HY and HO errors at all levels: {same} of {len(mols)}")

# ---- applicability survey (data/feasibility.csv, from feasibility_scan.py)
if os.path.exists("data/feasibility.csv"):
    import collections

    f = list(csv.DictReader(open("data/feasibility.csv", encoding="utf-8")))
    print("\n== applicability survey:", len({r['inchi'] for r in f}), "ATcT species with 3-5 heavy atoms")
    for s in ("ID", "HY", "HO"):
        c = collections.Counter(r["outcome"] for r in f if r["scheme"] == s)
        n = sum(c.values())
        print(f"  {s}: reaction {c['reaction']} ({100 * c['reaction'] / n:.0f}%), identity only {c['identity only']}, no reaction {c['no reaction']}, errors {c['error']}, timeouts {c['timeout']}")
    for g in ("N-containing", "halogen", "carbonyl", "ring", "nitrile", "amide", "aromatic"):
        sel = [r for r in f if r["scheme"] == "HO" and g in r["groups"].split(";")]
        c = collections.Counter(r["outcome"] for r in sel)
        print(f"  HO, {g:13s} n={len(sel):3d}: reaction {c['reaction']:3d} identity only {c['identity only']:3d} no reaction {c['no reaction']:3d}")
    t = [float(r["seconds"]) for r in f]
    print(f"  solve time: median {np.median(t) * 1000:.0f} ms, max {max(t):.2f} s")
