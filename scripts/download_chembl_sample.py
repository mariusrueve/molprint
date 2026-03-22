#!/usr/bin/env python3
"""Download a 10k SMILES sample from ChEMBL for corpus validation.

Usage:
    python scripts/download_chembl_sample.py [--size N] [--seed S]

Output: scripts/chembl_10k.smi  (SMILES TAB ChEMBL-ID, one per line)
"""

import argparse
import gzip
import os
import random
import sys
import urllib.request

CHEMBL_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    "chembl_36_chemreps.txt.gz"
)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTPUT = os.path.join(SCRIPT_DIR, "chembl_10k.smi")
DEFAULT_SIZE = 10_000


def download_sample(url: str, size: int, seed: int, output: str) -> None:
    print(f"Downloading ChEMBL SMILES from {url}")
    print("(This may take several minutes for the full file.)")

    random.seed(seed)
    reservoir: list[tuple[str, str]] = []
    total_seen = 0

    with urllib.request.urlopen(url) as resp:
        with gzip.open(resp, "rt", encoding="utf-8") as f:
            header = f.readline()  # skip header row
            print(f"Header: {header.strip()}")

            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                chembl_id = parts[0]
                smiles = parts[1]
                if not smiles or smiles == "None":
                    continue

                total_seen += 1
                if len(reservoir) < size:
                    reservoir.append((chembl_id, smiles))
                else:
                    j = random.randint(0, total_seen - 1)
                    if j < size:
                        reservoir[j] = (chembl_id, smiles)

    print(f"Scanned {total_seen} records, sampled {len(reservoir)}")

    with open(output, "w") as f:
        for chembl_id, smiles in reservoir:
            f.write(f"{smiles}\t{chembl_id}\n")

    print(f"Saved {len(reservoir)} SMILES to {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=DEFAULT_SIZE)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    download_sample(args.url if hasattr(args, "url") else CHEMBL_URL,
                    args.size, args.seed, args.output)


if __name__ == "__main__":
    main()
