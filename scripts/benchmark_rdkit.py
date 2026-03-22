#!/usr/bin/env python3
"""Measure RDKit fingerprint throughput for comparison against molprint.

Output: rdkit_bench.txt  (tab-separated: metric, molecules/sec)

Usage:
    python scripts/benchmark_rdkit.py [--molecules N] [--output FILE]

Requires: RDKit
"""

import argparse
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent

# Representative drug-like SMILES (diverse, realistic sizes)
BENCH_SMILES = [
    "CC(=O)Oc1ccccc1C(=O)O",                        # aspirin
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",          # testosterone
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",                 # caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                   # ibuprofen
    "c1ccc2c(c1)ccc3cccc4ccc5ccccc5c4c32",           # pyrene
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",       # glucose
    "CC(=O)NC1=CC=C(C=C1)O",                         # paracetamol
    "OC(=O)c1ccccc1O",                               # salicylic acid
    "CCOC(=O)c1ccc(cc1)N",                           # benzocaine
    "CN(C)CCC=c1[nH]c2ccccc2c1=C",                  # imipramine
    "Clc1ccc(cc1)C(c1ccc(Cl)cc1)C(Cl)(Cl)Cl",       # DDT
    "CC1=CC2=C(C=C1C)N(C=N2)CC(O)=O",               # similar to riboflavin fragment
    "O=C(O)c1ccncc1",                                # nicotinic acid
    "c1ccc(cc1)N",                                   # aniline
    "c1ccc2[nH]ccс2c1",                             # indole-like (may fail, skip)
    "O=C1NC(=O)NC=1",                                # uracil
]

# Filter to only parseable SMILES
def get_valid_smiles(Chem, smiles_list):
    valid = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            valid.append((smi, mol))
    return valid


def benchmark_morgan(Chem, AllChem, mols, n_repeat: int, radius: int, nbits: int) -> float:
    """Returns molecules/second."""
    t0 = time.perf_counter()
    for _ in range(n_repeat):
        for _, mol in mols:
            AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    elapsed = time.perf_counter() - t0
    total = n_repeat * len(mols)
    return total / elapsed


def benchmark_maccs(MACCSkeys, mols, n_repeat: int) -> float:
    """Returns molecules/second."""
    t0 = time.perf_counter()
    for _ in range(n_repeat):
        for _, mol in mols:
            MACCSkeys.GenMACCSKeys(mol)
    elapsed = time.perf_counter() - t0
    total = n_repeat * len(mols)
    return total / elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat", type=int, default=1000,
                        help="Number of repeat passes over the molecule set")
    parser.add_argument("--output", default=str(SCRIPT_DIR / "rdkit_bench.txt"))
    args = parser.parse_args()

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, MACCSkeys
    except ImportError:
        sys.exit("RDKit not found. Install with: conda install -c conda-forge rdkit")

    mols = get_valid_smiles(Chem, BENCH_SMILES)
    if not mols:
        sys.exit("No valid SMILES parsed.")
    print(f"Benchmarking RDKit on {len(mols)} molecules × {args.repeat} repeats")

    results: list[tuple[str, float]] = []

    # Morgan r=2, 2048 bits
    tps = benchmark_morgan(Chem, AllChem, mols, args.repeat, radius=2, nbits=2048)
    print(f"  Morgan(r=2, 2048): {tps:,.0f} mol/s")
    results.append(("morgan_r2_2048", tps))

    # Morgan r=3, 2048 bits
    tps = benchmark_morgan(Chem, AllChem, mols, args.repeat, radius=3, nbits=2048)
    print(f"  Morgan(r=3, 2048): {tps:,.0f} mol/s")
    results.append(("morgan_r3_2048", tps))

    # MACCS 166
    tps = benchmark_maccs(MACCSkeys, mols, args.repeat)
    print(f"  MACCS-166:         {tps:,.0f} mol/s")
    results.append(("maccs_166", tps))

    output = Path(args.output)
    with open(output, "w") as f:
        f.write("metric\tmol_per_sec\n")
        for metric, tps in results:
            f.write(f"{metric}\t{tps:.1f}\n")
    print(f"\nResults written to {output}")


if __name__ == "__main__":
    main()
