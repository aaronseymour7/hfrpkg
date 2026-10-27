"""Publication figures for the manuscript.  Run from the paper/ directory:

    python make_figures.py          # needs data/benchmark_errors.csv and data/scaling.csv

Writes figures/*.png (600 dpi) and figures/*.pdf (vector).
"""
import csv
import io
import math
import os
import random
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 7.5, "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3, "ytick.major.size": 3, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none", "figure.dpi": 150,
})
# Okabe-Ito colour-blind-safe palette
COL = {"ID": "#E69F00", "HY": "#56B4E9", "HO": "#009E73"}
NAME = {"ID": "Isodesmic", "HY": "Hypohomodesmotic", "HO": "Homodesmotic"}
SCHEMES = ["ID", "HY", "HO"]
LEVELS = ["B3LYP/6-31+G(d)", "M06-2X/6-31+G(d)", "wB97X-D/def2-TZVP"]
LEVEL_LABEL = {"B3LYP/6-31+G(d)": "B3LYP/6-31+G(d)", "M06-2X/6-31+G(d)": "M06-2X/6-31+G(d)", "wB97X-D/def2-TZVP": "ωB97X-D/def2-TZVP"}
OUT = "figures"
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------- data
rows = list(csv.DictReader(open("data/benchmark_errors.csv", encoding="utf-8")))
mol_order = []
for r in rows:
    if r["inchi"] not in mol_order:
        mol_order.append(r["inchi"])
CLASS = {r["inchi"]: r["class"] for r in rows}
SMILES = {r["inchi"]: r["smiles"] for r in rows}
# err[stage][level][scheme] -> array over molecules (in mol_order)
err = {st: {lv: {s: np.zeros(len(mol_order)) for s in SCHEMES} for lv in LEVELS} for st in ("opt", "sp")}
for r in rows:
    i = mol_order.index(r["inchi"])
    for s in SCHEMES:
        err[r["stage"]][r["geometry_level"]][s][i] = float(r[s])
is_str = np.array([CLASS[m] == "strained" for m in mol_order])
SUBSETS = {"unstrained": ~is_str, "strained": is_str, "all": np.ones(len(mol_order), bool)}


def pooled(stage, scheme, subset="all"):
    m = SUBSETS[subset]
    return np.concatenate([err[stage][lv][scheme][m] for lv in LEVELS])


rmse = lambda v: float(np.sqrt(np.mean(np.square(v))))
mae = lambda v: float(np.mean(np.abs(v)))


def boot_rmse(stage, scheme, subset, n=4000, seed=11):
    """95% CI of the pooled RMSE, resampling molecules (all three levels of a molecule kept together)."""
    rng = np.random.default_rng(seed)
    m = np.where(SUBSETS[subset])[0]
    stack = np.stack([err[stage][lv][scheme] for lv in LEVELS])  # 3 x mols
    vals = []
    for _ in range(n):
        idx = rng.choice(m, size=len(m), replace=True)
        vals.append(math.sqrt(np.mean(np.square(stack[:, idx]))))
    return np.percentile(vals, [2.5, 97.5])


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight", pad_inches=0.04, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    print("wrote", name)


def panel(ax, letter, dx=-0.12, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")


def mol_png(smi, w=420, h=420, label_ends=True):
    m = Chem.MolFromSmiles(smi)
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    o = d.drawOptions()
    o.clearBackground = False
    o.bondLineWidth = 2.4
    o.padding = 0.10
    o.useBWAtomPalette()
    o.fixedBondLength = 40
    if label_ends:  # CH3 / CH2 labels make tiny species (ethane, ethylene) readable
        for a in m.GetAtoms():
            if a.GetSymbol() == "C" and a.GetDegree() == 1:
                n = a.GetTotalNumHs()
                o.atomLabels[a.GetIdx()] = "CH" + (str(n) if n > 1 else "")
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    return np.array(Image.open(io.BytesIO(d.GetDrawingText())).convert("RGBA"))


def crop(img, pad=6):
    """Trim transparent margins so structures can be placed at a predictable size."""
    a = img[..., 3] > 0
    ys, xs = np.where(a)
    return img[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]


def common_scale(smiles, box_w, box_h, **kw):
    """Inches-per-pixel that fits every structure in the box, so all share one bond length."""
    out = []
    for smi in smiles:
        ih, iw = crop(mol_png(smi, **kw)).shape[:2]
        out.append(min(box_w / iw, box_h / ih))
    return min(out)


def place_mol(fig, smi, cx, cy, box_w, box_h, scale=None, **kw):
    """Draw a structure centred at (cx, cy) inches with correct aspect; `scale` (in/px) is shared across a figure."""
    img = crop(mol_png(smi, **kw))
    ih, iw = img.shape[:2]
    scale = scale or min(box_w / iw, box_h / ih)
    w, h = iw * scale, ih * scale
    W, H = fig.get_size_inches()
    ax = fig.add_axes([(cx - w / 2) / W, (cy - h / 2) / H, w / W, h / H])
    ax.imshow(img, aspect="equal", interpolation="lanczos")
    ax.axis("off")
    return ax


# ============================================================ wildcard fragments drawn as Kekule structures
def frag_png(pattern, syntax, w=520, h=300):
    """A wildcard pattern drawn as an explicit Kekule structure: every atom is a letter, bonds are lines.
    X = any atom, Y = any non-carbon atom ([!#6]); hydrogens are not shown (they are added to complete valence)."""
    from rdkit.Chem import rdDepictor

    m = Chem.MolFromSmarts(pattern) if syntax == "SMARTS" else Chem.MolFromSmiles(pattern)
    if syntax != "SMARTS":
        Chem.RemoveStereochemistry(m)
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    rdDepictor.SetPreferCoordGen(False)
    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    o = d.drawOptions()
    o.clearBackground = False
    o.bondLineWidth = 2.4
    o.padding = 0.12
    o.fixedBondLength = 42
    o.useBWAtomPalette()
    for a in m.GetAtoms():
        sm = a.GetSmarts() if syntax == "SMARTS" else ""
        if a.GetAtomicNum() == 0:
            o.atomLabels[a.GetIdx()] = "Y" if "!#6" in sm else "X"
        else:  # e.g. the H atoms of H2: label every atom explicitly
            o.atomLabels[a.GetIdx()] = a.GetSymbol()
    d.DrawMolecule(m)
    d.FinishDrawing()
    return np.array(Image.open(io.BytesIO(d.GetDrawingText())).convert("RGBA"))


def fig_wildcards(keep, outname, isodesmic_note):
    rows = list(csv.DictReader(open("data/wildcard_libraries.csv", encoding="utf-8")))
    levels = [l for l in [(1, "Isogyric"), (2, "Isodesmic (candidates built directly from the bond types of the target)"), (3, "Hypohomodesmotic"), (4, "Homodesmotic")] if l[0] in keep]
    rows = [r for r in rows if int(r["level"]) in keep]
    ncol, cw, ch, W = 6, 1.07, 0.68, 6.6
    blocks = []  # (level title, [(role, [rows...])])
    for lv, name in levels:
        sel = [r for r in rows if int(r["level"]) == lv]
        blocks.append((lv, name, [("Co-reactants", [r for r in sel if r["role"] == "reactants"]),
                                  ("Candidate products", [r for r in sel if r["role"] == "products"])]))
    H = 0.62 if isodesmic_note else 0.45  # room for the legend
    for _, _, groups in blocks:
        H += 0.32
        for _, items in groups:
            H += 0.22 + math.ceil(len(items) / ncol) * ch
        H += 0.12
    imgs = {(r["level"], r["role"], r["pattern"]): frag_png(r["pattern"], r["syntax"]) for r in rows}
    sc = min(min(0.9 * cw / crop(i).shape[1], 0.50 * ch / crop(i).shape[0]) for i in imgs.values())
    fig = plt.figure(figsize=(W, H))
    y = H - 0.1
    for lv, name, groups in blocks:  # noqa: E501
        fig.text(0.05 / W, (y - 0.14) / H, f"Level {lv}: {name}", fontsize=9.5, fontweight="bold", va="center")
        fig.add_artist(plt.Line2D([0.05 / W, (W - 0.05) / W], [(y - 0.27) / H] * 2, color="#999999", lw=0.6, transform=fig.transFigure))
        y -= 0.32
        for gname, items in groups:
            fig.text(0.12 / W, (y - 0.09) / H, gname, fontsize=7.5, style="italic", color="#444444", va="center")
            y -= 0.22
            for i, r in enumerate(items):
                col, row = i % ncol, i // ncol
                cx = 0.12 + (col + 0.5) * cw
                cy = y - row * ch - 0.30 * ch
                img = crop(imgs[(r["level"], r["role"], r["pattern"])])
                iw, ih = img.shape[1] * sc, img.shape[0] * sc
                ax = fig.add_axes([(cx - iw / 2) / W, (cy - ih / 2) / H, iw / W, ih / H])
                ax.imshow(img, aspect="equal", interpolation="lanczos")
                ax.axis("off")
                fig.text(cx / W, (y - row * ch - 0.80 * ch) / H, r["pattern"].replace("*", "\u2217"), fontsize=5.4, family="DejaVu Sans Mono",
                         ha="center", va="center", color="#555555")
            y -= math.ceil(len(items) / ncol) * ch
        y -= 0.12
    fig.text(0.05 / W, 0.26 / H, "X: any atom;  Y: any non-carbon atom;  hydrogens are added to complete the valence.", fontsize=6.4, color="#444444", va="bottom")
    if isodesmic_note:
        fig.text(0.05 / W, 0.10 / H, "Isodesmic: X stands for the elements of the target; one X–X', X=X' or X≡X' species per bond type present, the hydride of each element, and H$_2$.",
                 fontsize=6.4, color="#444444", va="bottom")
    save(fig, outname)


# ============================================================ Figure 1: the three schemes on one molecule
def fig1():
    from hfrpkg.core import Homodesmotic, Hypohomodesmotic, Isodesmic

    target = "CCC(C)=O"  # butan-2-one
    fns = {"ID": Isodesmic, "HY": Hypohomodesmotic, "HO": Homodesmotic}
    rxn = {}
    inp = Chem.MolToInchi(Chem.AddHs(Chem.MolFromSmiles(target)))
    for k, f in fns.items():
        rhs, lhs, st = f(Chem.AddHs(Chem.MolFromSmiles(target)), None, None, None, None)
        assert st == "Optimal", k
        L = [(c, Chem.MolToSmiles(Chem.RemoveHs(m)), Chem.MolToInchi(m) == inp) for m, c in lhs]
        R = [(c, Chem.MolToSmiles(Chem.RemoveHs(m)), False) for m, c in rhs]
        L.sort(key=lambda t: not t[2])  # target first
        rxn[k] = (L, R)
    conserve = {
        "ID": "conserves bond types\n(atom pair and bond order)",
        "HY": "conserves heavy atoms by\nhybridization and H count",
        "HO": "conserves bond hybridization\ntypes and H environments",
    }
    W, rowH, hdrW = 6.6, 1.35, 1.55
    fig = plt.figure(figsize=(W, rowH * 3 + 0.1))
    sc = common_scale([t[1] for k in SCHEMES for side in rxn[k] for t in side], 0.96, 0.98)
    H = fig.get_size_inches()[1]
    cell, opw, coefw = 1.06, 0.26, 0.0
    for r, k in enumerate(SCHEMES):
        y0 = H - (r + 1) * rowH - 0.05 + rowH / 2 - 0.02  # row centre
        fig.add_artist(plt.Rectangle((0.0, (y0 - 0.5) / H), 0.06 / W, 1.0 / H, transform=fig.transFigure, color=COL[k], lw=0))
        fig.text(0.16 / W, (y0 + 0.28) / H, NAME[k], fontsize=9.5, fontweight="bold", va="center")
        fig.text(0.16 / W, (y0 - 0.20) / H, conserve[k], fontsize=7, va="center", color="#444444", linespacing=1.3)
        L, R = rxn[k]
        items = []
        for j, (c, s_, tgt) in enumerate(L):
            items += [("mol", c, s_, tgt)] + ([("op", "+")] if j < len(L) - 1 else [])
        items.append(("op", "→"))
        for j, (c, s_, tgt) in enumerate(R):
            items += [("mol", c, s_, False)] + ([("op", "+")] if j < len(R) - 1 else [])
        total = sum(cell if it[0] == "mol" else opw for it in items)
        x = hdrW + (W - hdrW - total) / 2
        for it in items:
            if it[0] == "op":
                fig.text((x + opw / 2) / W, y0 / H, it[1], ha="center", va="center", fontsize=13)
                x += opw
            else:
                _, c, s_, tgt = it
                place_mol(fig, s_, x + cell / 2, y0 - 0.02, cell - 0.06, 0.98, scale=sc)
                if c > 1:
                    fig.text((x + 0.02) / W, (y0 + 0.40) / H, f"{c} ×", fontsize=9, fontweight="bold", va="center")
                if tgt:
                    fig.add_artist(plt.Rectangle((x / W, (y0 - 0.5) / H), cell / W, 1.0 / H, transform=fig.transFigure,
                                                 fill=False, ec=COL[k], lw=1.1, ls=(0, (3, 2))))
                x += cell
    save(fig, "fig1_schemes")
    return rxn


# ============================================================ Figure 2: error distributions
def kde_profile(v, factor=0.35, cut=2.0, n=400):
    """Gaussian KDE (bandwidth = factor x sample SD, as in matplotlib's violinplot) evaluated `cut`
    bandwidths beyond the data, so the violin tapers to a point instead of ending in a flat cap at the
    extreme observations.  Returns (y, density scaled to a maximum of 1)."""
    v = np.asarray(v, float)
    h = factor * v.std(ddof=1)
    y = np.linspace(v.min() - cut * h, v.max() + cut * h, n)
    d = np.exp(-0.5 * ((y[:, None] - v[None, :]) / h) ** 2).sum(axis=1)
    return y, d / d.max()


def dist_panel(ax, stage, subset, ylim):
    positions = np.arange(3)
    for p, s in zip(positions, SCHEMES):
        v = pooled(stage, s, subset)
        ky, kd = kde_profile(v)
        ax.fill_betweenx(ky, p - 0.425 * kd, p + 0.425 * kd, facecolor=COL[s], alpha=0.28, edgecolor="none", lw=0)
        ax.plot(np.r_[p - 0.425 * kd, (p + 0.425 * kd)[::-1]], np.r_[ky, ky[::-1]], color=COL[s], alpha=0.55, lw=0.8)
        rng = np.random.default_rng(3 + p)
        ax.scatter(p + rng.uniform(-0.2, 0.2, len(v)), v, s=7, color=COL[s], alpha=0.65, lw=0, zorder=3)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.add_patch(plt.Rectangle((p - 0.07, q1), 0.14, q3 - q1, fc="white", ec="black", lw=0.7, zorder=4))
        ax.plot([p - 0.07, p + 0.07], [med, med], color="black", lw=1.4, zorder=5)
    ax.axhline(0, color="black", lw=0.7, ls=(0, (4, 3)), zorder=1)
    ax.axhspan(-1, 1, color="#999999", alpha=0.10, lw=0, zorder=0)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{NAME[s]}\n({s})" if False else s for s in SCHEMES])
    ax.set_ylim(*ylim)
    ax.set_xlim(-0.6, 2.6)
    for p, s in zip(positions, SCHEMES):
        v = pooled(stage, s, subset)
        ax.text(p, ylim[1] - 0.05, f"RMSE {rmse(v):.2f}\nMAE {mae(v):.2f}", ha="center", va="top", fontsize=7.2,
                bbox=dict(fc="white", ec="none", pad=1.2, alpha=0.85))


def fig2(stage="sp", name="fig2_distributions", ylim=None):
    if ylim is None:  # derived from the data so that no point is ever clipped
        ext = [kde_profile(pooled(stage, s, sub))[0] for s in SCHEMES for sub in ("unstrained", "strained")]
        ylim = (np.floor(min(e.min() for e in ext) - 0.3), np.ceil(max(e.max() for e in ext) + 2.4))
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.5), sharey=True, gridspec_kw=dict(wspace=0.06, width_ratios=[1, 1]))
    for ax, subset, title in zip(axes, ("unstrained", "strained"), ("Unstrained (acyclic, {} molecules)", "Strained (3–5-membered rings, {} molecules)")):
        dist_panel(ax, stage, subset, ylim)
        n = int(SUBSETS[subset].sum()) * 3
        title = title.format(n // 3)
        ax.set_title(f"{title}\n{n} reactions", fontsize=8.2, pad=6)
    axes[0].set_ylabel("Error vs. ATcT (kcal mol$^{-1}$)")
    axes[1].spines["left"].set_visible(False)
    axes[1].tick_params(axis="y", length=0)
    panel(axes[0], "a", -0.16)
    panel(axes[1], "b", -0.04)
    save(fig, name)


# ============================================================ Figure 3: RMSE before/after SPR with bootstrap CIs
def fig3():
    fig, axes = plt.subplots(1, 3, figsize=(6.6, 3.0), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, subset, title in zip(axes, ("unstrained", "strained", "all"), ("Unstrained", "Strained", "All molecules")):
        for j, s in enumerate(SCHEMES):
            for k, (stage, hatch, alpha) in enumerate((("opt", "////", 0.45), ("sp", None, 1.0))):
                v = pooled(stage, s, subset)
                lo, hi = boot_rmse(stage, s, subset)
                x = j + (k - 0.5) * 0.38
                ax.bar(x, rmse(v), width=0.36, color=COL[s], alpha=alpha, hatch=hatch, edgecolor="white" if hatch is None else COL[s], lw=0)
                ax.errorbar(x, rmse(v), yerr=[[rmse(v) - lo], [hi - rmse(v)]], color="black", lw=0.8, capsize=2, capthick=0.8)
                ax.text(x, hi + 0.08, f"{rmse(v):.2f}", ha="center", va="bottom", fontsize=6.6)
        ax.set_xticks(range(3)); ax.set_xticklabels(SCHEMES)
        ax.set_title(title, fontsize=8.5, pad=8)
        ax.set_ylim(0, 5.2)
        ax.set_xlim(-0.6, 2.6)
    axes[0].set_ylabel("RMSE vs. ATcT (kcal mol$^{-1}$)")
    for a_ in axes[1:]:
        a_.spines["left"].set_visible(False); a_.tick_params(axis="y", length=0)
    for a_, l in zip(axes, "abc"):
        panel(a_, l, -0.02 if a_ is not axes[0] else -0.17, 1.06)
    handles = [Patch(fc="#888888", alpha=0.45, hatch="////", ec="#888888", label="Before SPR (lower-level energies)"),
               Patch(fc="#888888", label="After SPR (ωB97X-V/def2-TZVP)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.08))
    save(fig, "fig3_rmse")


# ============================================================ Figure 4: geometry independence
def fig4():
    fig = plt.figure(figsize=(6.6, 2.95))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.12], wspace=0.55, left=0.09, right=0.99, top=0.88, bottom=0.2)
    lv_a, lv_b = LEVELS[0], LEVELS[2]
    for i, (stage, title) in enumerate((("opt", "Before SPR"), ("sp", "After SPR"))):
        ax = fig.add_subplot(gs[i])
        x = err[stage][lv_a]["HO"]; y = err[stage][lv_b]["HO"]
        lim = float(np.ceil(max(np.abs(x).max(), np.abs(y).max()) + 0.5))
        ax.plot([-lim, lim], [-lim, lim], color="black", lw=0.6, ls=(0, (4, 3)), zorder=1)
        ax.scatter(x[~is_str], y[~is_str], s=14, color=COL["HO"], lw=0, alpha=0.8, label="unstrained", zorder=3)
        ax.scatter(x[is_str], y[is_str], s=16, color="white", edgecolor=COL["HO"], lw=0.9, label="strained", zorder=3)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
        ax.set_title(title, fontsize=8.5, pad=6)
        ax.set_xlabel("Error, B3LYP geometry")
        ax.set_ylabel("Error, ωB97X-D geometry")
        if i == 0:
            ax.legend(frameon=False, loc="upper left", handletextpad=0.2, borderpad=0, labelspacing=0.2, fontsize=7)
        ax.text(0.97, 0.04, f"max |Δ| = {np.max(np.abs(x - y)):.2f}", transform=ax.transAxes, ha="right", fontsize=6.8)
        panel(ax, "ab"[i], -0.30, 1.05)
    ax = fig.add_subplot(gs[2])
    for s in SCHEMES:
        for stage, ls in (("opt", (0, (3, 2))), ("sp", "-")):
            stack = np.array([err[stage][lv][s] for lv in LEVELS])
            xs = np.sort(stack.max(0) - stack.min(0))
            ax.step(np.r_[0, xs], np.r_[0, np.arange(1, len(xs) + 1) / len(xs)], where="post", color=COL[s], ls=ls, lw=1.5)
    ax.set_xlim(0, 4)
    ax.set_xlabel("Spread across 3 geometries\n(kcal mol$^{-1}$)")
    ax.set_ylabel("Cumulative fraction of molecules")
    ax.set_ylim(0, 1.02)
    h = [Line2D([], [], color=COL[s], lw=1.6, label=s) for s in SCHEMES]
    h += [Line2D([], [], color="black", lw=1.2, ls=(0, (3, 2)), label="before SPR"), Line2D([], [], color="black", lw=1.2, label="after SPR")]
    ax.legend(handles=h, frameon=False, loc="lower right", handlelength=1.8, labelspacing=0.25, fontsize=7)
    panel(ax, "c", -0.26, 1.05)
    save(fig, "fig4_geometry")


# ============================================================ Figure 5: strained molecules
def fig5():
    idx = np.where(is_str)[0]
    mean = {s: np.mean([err["sp"][lv][s][idx] for lv in LEVELS], axis=0) for s in SCHEMES}
    order = np.argsort(mean["HO"])[::-1]  # largest positive error at top
    n = len(order)
    W, H = 6.6, 0.38 * n + 0.9
    fig = plt.figure(figsize=(W, H))
    left, right, bottom, top = 1.3, 0.15, 0.62, 0.15
    sc = common_scale([SMILES[mol_order[idx[o]]] for o in order], 1.1, 0.34, label_ends=False)
    ax = fig.add_axes([left / W, bottom / H, (W - left - right) / W, (H - bottom - top) / H])
    y = np.arange(n)[::-1]
    ax.axvspan(-1, 1, color="#999999", alpha=0.10, lw=0)
    ax.axvline(0, color="black", lw=0.7, ls=(0, (4, 3)))
    for yi, o in zip(y, order):
        vals = [mean[s][o] for s in SCHEMES]
        ax.plot([min(vals), max(vals)], [yi, yi], color="#BBBBBB", lw=1.0, zorder=1)
    for s, mk in zip(SCHEMES, ("s", "^", "o")):
        ax.scatter([mean[s][o] for o in order], y, marker=mk, s=30, color=COL[s], edgecolor="black", lw=0.4, zorder=3, label=NAME[s])
    ax.set_yticks([]); ax.set_ylim(-0.7, n - 0.3); ax.set_xlim(-5.5, 8.5)
    ax.set_xlabel("Error vs. ATcT after SPR (kcal mol$^{-1}$; mean of three geometries)")
    ax.spines["left"].set_visible(False)
    ax.legend(frameon=False, loc="lower right", handletextpad=0.2)
    ylo, yhi = ax.get_ylim()
    for yi, o in zip(y, order):
        # structure thumbnail in the left margin, aligned with its row
        cy = bottom + (yi - ylo) / (yhi - ylo) * (H - bottom - top)
        place_mol(fig, SMILES[mol_order[idx[o]]], left / 2 - 0.05, cy, 1.1, 0.34, scale=sc, label_ends=False)
    save(fig, "fig5_strained")


# ============================================================ Figure 6: scalability
def fig6():
    data = list(csv.DictReader(open("data/scaling.csv")))
    W, H = 6.6, 3.85
    fig = plt.figure(figsize=(W, H))
    axa = fig.add_axes([0.085, 0.46, 0.40, 0.49]); axb = fig.add_axes([0.585, 0.46, 0.395, 0.49])
    mk = {"polyene": ("o", "#0072B2"), "PAH": ("s", "#D55E00"), "fullerene": ("*", "#CC79A7")}
    lab = {"polyene": "Linear polyenes (C$_6$–C$_{38}$)", "PAH": "Benzene to coronene", "fullerene": "C$_{60}$"}
    for fam in ("polyene", "PAH", "fullerene"):
        d = [r for r in data if r["family"] == fam]
        x = [int(r["heavy_atoms"]) for r in d]
        m, c = mk[fam]
        ms = 130 if fam == "fullerene" else 28
        axa.scatter(x, [int(r["distinct_refs"]) for r in d], marker=m, s=ms, color=c, edgecolor="black", lw=0.4, label=lab[fam], zorder=3)
        axb.scatter(x, [int(r["total_ref_molecules"]) for r in d], marker=m, s=ms, color=c, edgecolor="black", lw=0.4, zorder=3)
    axa.set_ylim(0, 8); axa.set_yticks([0, 2, 4, 6, 8])
    axa.set_xlabel("Heavy atoms in target"); axa.set_ylabel("Distinct reference species")
    axa.legend(frameon=False, loc="upper left", fontsize=7, handletextpad=0.1)
    axb.set_xlabel("Heavy atoms in target"); axb.set_ylabel("Reference molecules (Σ coefficients)")
    mx = max(float(r["solve_s"]) for r in data)
    axb.text(0.04, 0.93, f"ILP solve time ≤ {mx:.2f} s\nfor all {len(data)} targets", transform=axb.transAxes, fontsize=7, va="top")
    panel(axa, "a", -0.22, 1.04); panel(axb, "b", -0.22, 1.04)
    y0 = 0.62
    fig.text(0.012, (y0 + 0.45) / H, "c", fontsize=11, fontweight="bold", va="center")
    smis = ["CC=C", "C=C", "C=CC(=C)C"]
    sc = common_scale(smis, 0.85, 0.62)
    x = 0.45
    fig.text((x + 0.3) / W, y0 / H, "C$_{60}$", ha="center", va="center", fontsize=13, fontweight="bold"); x += 0.75
    for sym, smi, cnt in (("+", smis[0], "60 ×"), ("+", smis[1], "30 ×"), ("→", smis[2], "60 ×")):
        fig.text((x + 0.13) / W, y0 / H, sym, ha="center", va="center", fontsize=13); x += 0.28
        fig.text(x / W, y0 / H, cnt, ha="left", va="center", fontsize=9.5, fontweight="bold"); x += 0.52
        place_mol(fig, smi, x + 0.45, y0, 0.9, 0.62, scale=sc); x += 0.95
    fig.text(0.5, 0.10 / H, "Reference species: propene, ethylene and isoprene (all with ATcT values); the set does not grow with fullerene size.",
             ha="center", va="center", fontsize=7.2, color="#333333")
    save(fig, "fig6_scalability")


if __name__ == "__main__":
    fig1()
    fig2("sp", "fig2_distributions")
    fig2("opt", "figS1_distributions_preSPR")
    fig3()
    fig4()
    fig5()
    fig6()
    fig_wildcards((1, 2, 3), "figS2_wildcard_libraries", True)
    fig_wildcards((4,), "figS3_wildcard_library_homodesmotic", False)
