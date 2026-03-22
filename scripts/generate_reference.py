#!/usr/bin/env python3
"""
Generate reference fingerprints using RDKit for validating molprint output.

Usage:
    python scripts/generate_reference.py

Outputs:
    scripts/reference_maccs.tsv   — SMILES, name, set bit indices (RDKit MACCS)
    scripts/reference_morgan.tsv  — SMILES, name, pairwise Tanimoto similarities

Requires: rdkit (pip install rdkit)
"""

import sys
from itertools import combinations

try:
    from rdkit import Chem
    from rdkit.Chem import MACCSkeys
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
except ImportError:
    print("ERROR: RDKit not found. Install with: pip install rdkit", file=sys.stderr)
    sys.exit(1)

MOLECULES = [
    # Simple / functional groups
    ("methane",          "C"),
    ("ethane",           "CC"),
    ("ethanol",          "CCO"),
    ("acetic_acid",      "CC(=O)O"),
    ("methylamine",      "CN"),
    ("methanethiol",     "CS"),
    ("chloromethane",    "CCl"),
    ("fluoromethane",    "CF"),
    ("bromomethane",     "CBr"),
    ("iodomethane",      "CI"),
    ("acetaldehyde",     "CC=O"),
    ("acetone",          "CC(=O)C"),
    ("dimethylether",    "COC"),
    ("dimethylsulfide",  "CSC"),
    ("methylamine",      "CN"),
    ("dimethylamine",    "CNC"),
    ("trimethylamine",   "CN(C)C"),
    ("urea",             "NC(=O)N"),
    ("acetamide",        "CC(=O)N"),
    ("nitromethane",     "C[N+](=O)[O-]"),
    # Rings
    ("benzene",          "c1ccccc1"),
    ("cyclohexane",      "C1CCCCC1"),
    ("cyclopentane",     "C1CCCC1"),
    ("toluene",          "Cc1ccccc1"),
    ("phenol",           "Oc1ccccc1"),
    ("aniline",          "Nc1ccccc1"),
    ("pyridine",         "c1ccncc1"),
    ("furan",            "c1ccoc1"),
    ("thiophene",        "c1ccsc1"),
    ("pyrrole",          "c1cc[nH]c1"),
    ("imidazole",        "c1cnc[nH]1"),
    ("indole",           "c1ccc2[nH]ccc2c1"),
    ("naphthalene",      "c1ccc2ccccc2c1"),
    # Drug-like
    ("aspirin",          "CC(=O)Oc1ccccc1C(=O)O"),
    ("caffeine",         "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    ("acetaminophen",    "CC(=O)Nc1ccc(O)cc1"),
    ("ibuprofen",        "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("lidocaine",        "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    ("ciprofloxacin",    "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O"),
]


def generate_maccs_reference():
    rows = []
    skipped = []
    for name, smi in MOLECULES:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            skipped.append((name, smi))
            continue
        fp = MACCSkeys.GenMACCSKeys(mol)
        bits = list(fp.GetOnBits())
        rows.append((smi, name, bits))

    if skipped:
        print(f"WARNING: {len(skipped)} molecules failed to parse:", file=sys.stderr)
        for name, smi in skipped:
            print(f"  {name}: {smi}", file=sys.stderr)

    out = "scripts/reference_maccs.tsv"
    with open(out, "w") as f:
        f.write("smiles\tname\tbits\n")
        for smi, name, bits in rows:
            f.write(f"{smi}\t{name}\t{','.join(map(str, bits))}\n")
    print(f"Written {len(rows)} MACCS records to {out}")
    return rows


def generate_morgan_reference(maccs_rows):
    """Compute pairwise Tanimoto for all molecule pairs using Morgan r=2 2048-bit."""
    fps = []
    for smi, name, _ in maccs_rows:
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        fps.append((name, smi, fp))

    out = "scripts/reference_morgan_sims.tsv"
    with open(out, "w") as f:
        f.write("name_a\tname_b\tsmiles_a\tsmiles_b\ttanimoto\n")
        for (na, sa, fa), (nb, sb, fb) in combinations(fps, 2):
            sim = DataStructs.TanimotoSimilarity(fa, fb)
            f.write(f"{na}\t{nb}\t{sa}\t{sb}\t{sim:.6f}\n")
    print(f"Written {len(list(combinations(fps, 2)))} Morgan similarity pairs to {out}")


if __name__ == "__main__":
    rows = generate_maccs_reference()
    generate_morgan_reference(rows)
    print("Done. Run 'cargo test -p molprint-fp validate' to compare against molprint.")
