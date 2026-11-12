"""Bulk-update hfrpkg's ATcT library (hfrpkg/data/ATcT_lib.txt) from the ATcT API.

Run from the hfrpkg repository root:

    python update_atct_lib.py --inspect          # download and show the fields of the first records
    python update_atct_lib.py                    # dry run: report what would be added / changed
    python update_atct_lib.py --write            # add species that are not in the library yet
    python update_atct_lib.py --write --update-existing   # ALSO overwrite values of species already there

Only gas-phase, neutral species are taken, with their 0 K enthalpy of formation (the same kind of entry
the library already holds). Condensed-phase (cr, l, aq) and ionic entries, and the 298.15 K values, are skipped.
By default existing entries are never changed, so the paper's benchmark numbers stay reproducible.
A copy of the old file is saved as ATcT_lib.txt.bak before anything is written.
The download is cached in atct_all.json; delete it to fetch again.
"""
import argparse
import json
import pathlib
import re
import shutil
import urllib.error
import urllib.request

URL = "https://atct.anl.gov/api/v1/all/"
ROOT = pathlib.Path(__file__).resolve().parent
LIB = ROOT / "hfrpkg" / "data" / "ATcT_lib.txt"
CACHE = ROOT / "atct_all.json"

# candidate key names (matched case-insensitively, ignoring spaces/underscores/punctuation)
KEYS = {
    "inchi": ["inchi", "inchistring", "standardinchi"],
    "smiles": ["smiles", "canonicalsmiles"],
    "hf0": ["dhf0k", "deltahf0k", "hf0k", "dfh0k", "hf0", "enthalpy0k", "deltafh0k", "dhf0"],
    "unc": ["uncertainty", "unc", "dhfuncertainty", "error", "hf0kuncertainty", "uncertainty0k"],
    "name": ["name", "speciesname", "species"],
    "formula": ["formula", "molecularformula"],
    "id": ["atctid", "id", "speciesid", "atct_id"],
    "phase": ["phase", "state", "physicalstate", "aggregation", "speciesphase"],
}
norm = lambda k: re.sub(r"[^a-z0-9]", "", k.lower())


def fetch():
    if not CACHE.exists():
        req = urllib.request.Request(URL, headers={"User-Agent": "hfrpkg-atct-updater", "Accept": "application/json"})
        data, url = [], URL
        while url:  # follow pagination if the API pages its results
            try:
                with urllib.request.urlopen(req if url == URL else urllib.request.Request(url, headers=req.headers), timeout=120) as r:
                    page = json.load(r)
            except urllib.error.HTTPError as e:
                raise SystemExit(
                    f"The ATcT server refused the download (HTTP {e.code}).\n"
                    f"Open {url} in your web browser and save it as atct_all.json in\n{ROOT}\n"
                    "then run this again; the script reads that file instead of downloading.")
            if isinstance(page, dict) and isinstance(page.get("results"), list):
                data += page["results"]; url = page.get("next")
            else:
                data = page; url = None
        CACHE.write_text(json.dumps(data), encoding="utf-8")
    data = json.loads(CACHE.read_text(encoding="utf-8"))
    if isinstance(data, dict):  # {"data": [...]} or {id: record}
        for k in ("data", "species", "results", "items"):
            if isinstance(data.get(k), list):
                return data[k]
        return [dict(v, id=k) if isinstance(v, dict) else v for k, v in data.items()]
    return data


def flat(rec, prefix=""):
    out = {}
    for k, v in rec.items():
        if isinstance(v, dict):
            out.update(flat(v, prefix + k + "_"))
        else:
            out[norm(prefix + k)] = v
    return out


def pick(rec, field):
    for k in KEYS[field]:
        if k in rec and rec[k] not in (None, ""):
            return rec[k]
    return None


def is_gas(r):
    ph = pick(r, "phase")
    if ph is not None:
        return str(ph).strip().lower().strip("()") in ("g", "gas", "gaseous", "gas phase", "ideal gas")
    tag = " ".join(str(pick(r, f) or "") for f in ("name", "formula")).lower()
    m = re.search(r"\((g|cr|cr,[^)]*|l|aq|s|liq|sol)\)", tag)
    return m is None or m.group(1) == "g"


def num(x):
    try:
        return float(str(x).replace("±", " ").split()[0])
    except (ValueError, IndexError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inspect", action="store_true"); ap.add_argument("--write", action="store_true")
    ap.add_argument("--update-existing", action="store_true")
    ap.add_argument("--hf0-key", help="field holding dHf(0 K) if not detected"); ap.add_argument("--unc-key")
    ap.add_argument("--phase-key", help="field holding the phase/state if not detected")
    ap.add_argument("--units", choices=["kJ", "kcal"], default="kJ", help="units of the API values (library is kJ/mol)")
    a = ap.parse_args()
    recs = [flat(r) for r in fetch() if isinstance(r, dict)]
    print(f"{len(recs)} records downloaded")
    if a.inspect:
        for r in recs[:3]:
            print(json.dumps(r, indent=1, default=str)[:2000]); print("-" * 60)
        return
    if a.hf0_key: KEYS["hf0"].insert(0, norm(a.hf0_key))
    if a.unc_key: KEYS["unc"].insert(0, norm(a.unc_key))
    if a.phase_key: KEYS["phase"].insert(0, norm(a.phase_key))
    if pick(recs[0], "phase") is None:
        print("Note: no phase field found; using the (g)/(cr)/(l)/(aq) tag in the name or formula instead.")
    if pick(recs[0], "hf0") is None:
        raise SystemExit("Could not find the 0 K enthalpy field. Run with --inspect and pass --hf0-key <name>.")

    from rdkit import Chem, RDLogger
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    RDLogger.DisableLog("rdApp.*")
    lines = LIB.read_text(encoding="utf-8").splitlines()
    have = {l.split("\t")[3]: i for i, l in enumerate(lines) if len(l.split("\t")) > 5}
    f = 4.184 if a.units == "kcal" else 1.0
    add, change, skip, notgas, ion = [], [], 0, 0, 0
    for r in recs:
        if not is_gas(r):
            notgas += 1; continue
        smi, inchi, hf0, unc = pick(r, "smiles"), pick(r, "inchi"), num(pick(r, "hf0")), num(pick(r, "unc"))
        if hf0 is None:
            skip += 1; continue
        mol = Chem.MolFromSmiles(smi) if smi else (Chem.MolFromInchi(inchi) if inchi else None)
        if mol is None:
            skip += 1; continue
        if Chem.GetFormalCharge(mol) != 0:
            ion += 1; continue
        key = Chem.MolToInchi(Chem.AddHs(mol))  # same key hfrpkg uses for look-ups
        form = pick(r, "formula") or "".join(e + (n or "1") for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", CalcMolFormula(mol)))
        row = "\t".join([str(pick(r, "id") or ""), str(pick(r, "name") or ""), str(form), key, Chem.MolToSmiles(mol),
                         f"{hf0 * f:.3f}", f"{(unc or 0) * f:.3f}"])
        if key not in have:
            add.append(row); have[key] = None
        elif have[key] is not None:
            old = lines[have[key]].split("\t")
            if abs(float(old[5]) - hf0 * f) > 1e-3:
                change.append((have[key], lines[have[key]], row))
    print(f"skipped: {notgas} not gas phase, {ion} ions, {skip} without a usable structure or 0 K value")
    print(f"new gas-phase neutral species: {len(add)}   existing with a different 0 K value: {len(change)}")
    for r in add[:15]: print("  + " + r)
    for _, o, n in change[:10]: print(f"  ~ {o.split(chr(9))[1]}: {o.split(chr(9))[5]} -> {n.split(chr(9))[5]} kJ/mol")
    if not a.write:
        print("Dry run; nothing written. Add --write to apply."); return
    shutil.copy(LIB, LIB.with_suffix(".txt.bak"))
    if a.update_existing:
        for i, _, n in change: lines[i] = n
    lines += add
    LIB.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {LIB} (+{len(add)} new" + (f", {len(change)} updated)" if a.update_existing else ", existing values unchanged)"))


if __name__ == "__main__":
    main()
