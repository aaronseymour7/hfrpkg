import glob
import os
import shutil

from hfrpkg.common import conformer_name, read_index, write_index


def numeric_key(path):
    head = os.path.basename(path).split(".")[0]
    return (0, int(head)) if head.isdigit() else (1, head)


def grab_unique_coms(software=None, output_dir="unique_files"):
    """Copy the input files of every distinct species (by InChI) across all *.mhfr into unique_files/.

    A species with several conformers contributes all of them (n-cNN files). Conformer sets are seeded
    per species, so every folder holds the same set and the first occurrence is representative.
    """
    os.makedirs(output_dir, exist_ok=True)
    seen = {}
    entries, counts = [], {}
    level_sw = software
    rule = None
    counter = 1

    for mhfr_dir in sorted(glob.glob("*.mhfr"), key=numeric_key):
        if not os.path.isfile(os.path.join(mhfr_dir, "index.txt")):
            continue
        index = read_index(mhfr_dir)
        level_sw = level_sw or index.software
        rule = rule or index.conformer_rule
        in_ext, _ = index.extensions
        for name, inchi, smiles in index.entries:
            if inchi in seen:
                continue
            n = index.conformers.get(name, 1)
            stems = [name] if n <= 1 else [conformer_name(name, k) for k in range(1, n + 1)]
            if not all(os.path.isfile(os.path.join(mhfr_dir, st + in_ext)) for st in stems):
                print(f"[WARNING] Missing input file(s) for {mhfr_dir}/{name}")
                continue
            new = str(counter)
            for k, st in enumerate(stems, start=1):
                dest = new if n <= 1 else conformer_name(new, k)
                shutil.copyfile(os.path.join(mhfr_dir, st + in_ext), os.path.join(output_dir, dest + in_ext))
            entries.append((new, inchi, smiles))
            if n > 1:
                counts[new] = n
            seen[inchi] = counter
            counter += 1

    write_index(output_dir, "Unique", level_sw, None, entries, conformers=counts, conformer_rule=rule)
    print(f"Done. {len(entries)} unique species ({sum(counts.values()) + len(entries) - len(counts)} input files) written to: {output_dir}/")
