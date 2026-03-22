#!/usr/bin/env python3
"""Validate molprint MACCS fingerprints against RDKit on a ChEMBL corpus.

Usage:
    python scripts/corpus_validate.py [smiles_file]

Requires:
    - RDKit  (conda install -c conda-forge rdkit)
    - molprint CLI built: cargo build --release

Exits 0 on success (>= 99% per-molecule MACCS accuracy), 1 on failure.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).parent
WORKSPACE = SCRIPT_DIR.parent
TARGET_ACCURACY = 99.0


def require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import MACCSkeys
        return Chem, MACCSkeys
    except ImportError:
        sys.exit(
            "RDKit not found.\n"
            "Install with:  conda install -c conda-forge rdkit"
        )


def find_molprint() -> Path:
    release = WORKSPACE / "target" / "release" / "molprint"
    if release.exists():
        return release
    found = shutil.which("molprint")
    if found:
        return Path(found)
    sys.exit(
        "molprint binary not found.\n"
        "Build with:  cargo build --release"
    )


def parse_fps_bits(hex_fp: str) -> set[int]:
    """Convert an FPS hex fingerprint string to a set of set bit indices."""
    bits: set[int] = set()
    for byte_idx in range(len(hex_fp) // 2):
        byte_val = int(hex_fp[byte_idx * 2 : byte_idx * 2 + 2], 16)
        for bit_offset in range(8):
            if byte_val & (1 << bit_offset):
                bits.add(byte_idx * 8 + bit_offset)
    return bits


def main() -> None:
    smi_file = Path(sys.argv[1]) if len(sys.argv) > 1 else SCRIPT_DIR / "chembl_10k.smi"
    if not smi_file.exists():
        sys.exit(
            f"SMILES file not found: {smi_file}\n"
            "Run scripts/download_chembl_sample.py first."
        )

    Chem, MACCSkeys = require_rdkit()
    molprint_bin = find_molprint()

    # ── Load SMILES ───────────────────────────────────────────────────────────
    smiles_list: list[str] = []
    ids_list: list[str] = []
    with open(smi_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            smiles_list.append(parts[0])
            ids_list.append(parts[1] if len(parts) > 1 else f"mol{len(ids_list)}")

    print(f"Loaded {len(smiles_list)} SMILES from {smi_file}")

    # ── RDKit MACCS ──────────────────────────────────────────────────────────
    rdkit_bits: dict[str, set[int]] = {}
    rdkit_failures = 0
    for smi, mid in zip(smiles_list, ids_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rdkit_failures += 1
            continue
        fp = MACCSkeys.GenMACCSKeys(mol)
        # RDKit MACCS bits are 1-indexed (bit 0 is unused, bits 1–165 are keys)
        rdkit_bits[mid] = {b for b in fp.GetOnBits() if 1 <= b <= 165}

    print(
        f"RDKit: {len(rdkit_bits)}/{len(smiles_list)} parsed "
        f"({rdkit_failures} failures)"
    )

    # ── molprint MACCS via CLI ────────────────────────────────────────────────
    with tempfile.NamedTemporaryFile(suffix=".fps", delete=False) as tmp:
        fps_tmp = Path(tmp.name)

    try:
        result = subprocess.run(
            [str(molprint_bin), "fp", "--fp-type", "maccs",
             "--input", str(smi_file), "--output", str(fps_tmp)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            print(f"molprint error:\n{result.stderr}", file=sys.stderr)
            sys.exit(1)

        molprint_bits: dict[str, set[int]] = {}
        with open(fps_tmp) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t", 1)
                if len(parts) < 2:
                    continue
                hex_fp, mol_id = parts[0], parts[1]
                # Molprint MACCS: bit N = MACCS key N (0-indexed); compare 1–165
                all_bits = parse_fps_bits(hex_fp)
                molprint_bits[mol_id] = {b for b in all_bits if 1 <= b <= 165}
    finally:
        fps_tmp.unlink(missing_ok=True)

    print(f"molprint: {len(molprint_bits)} fingerprints computed")

    # ── Compare ───────────────────────────────────────────────────────────────
    common = set(rdkit_bits) & set(molprint_bits)
    if not common:
        sys.exit("No common molecule IDs to compare — check ID formats.")

    mismatches = 0
    bit_errors: dict[int, list[int, int]] = {}  # bit → [fp, fn]
    for mid in sorted(common):
        rb = rdkit_bits[mid]
        mb = molprint_bits[mid]
        if rb != mb:
            mismatches += 1
            for b in mb - rb:
                bit_errors.setdefault(b, [0, 0])[0] += 1  # false positive
            for b in rb - mb:
                bit_errors.setdefault(b, [0, 0])[1] += 1  # false negative

    correct = len(common) - mismatches
    accuracy = 100.0 * correct / len(common)

    print(f"\nResults ({len(common)} molecules compared):")
    print(f"  Exact MACCS matches : {correct}  ({accuracy:.2f}%)")
    print(f"  Mismatches          : {mismatches}")

    if bit_errors:
        print("\nTop bit errors (bit: false-positives / false-negatives):")
        top = sorted(bit_errors.items(), key=lambda x: -(x[1][0] + x[1][1]))[:20]
        for bit, (fp, fn) in top:
            print(f"  key {bit:3d}: {fp} FP, {fn} FN")

    if accuracy >= TARGET_ACCURACY:
        print(f"\n✓ PASS: accuracy {accuracy:.2f}% >= {TARGET_ACCURACY}%")
    else:
        print(f"\n✗ FAIL: accuracy {accuracy:.2f}% < {TARGET_ACCURACY}%", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
