"""Smoke tests for the molprint Python bindings.

Run after `maturin develop`:
    python python/test_bindings.py
"""

import sys


def check(condition: bool, msg: str) -> None:
    if not condition:
        print(f"FAIL: {msg}", file=sys.stderr)
        sys.exit(1)
    print(f"  ok: {msg}")


def main() -> None:
    try:
        import molprint
    except ImportError as e:
        sys.exit(f"Cannot import molprint: {e}\nRun: maturin develop")

    print("molprint Python bindings smoke test")
    print("=" * 40)

    # ── Mol.from_smiles ────────────────────────────────────────────────────
    mol = molprint.Mol.from_smiles("c1ccccc1")
    check(mol is not None, "Mol.from_smiles('c1ccccc1') succeeds")

    try:
        molprint.Mol.from_smiles("not_valid!!!")
        check(False, "invalid SMILES should raise ValueError")
    except ValueError:
        check(True, "invalid SMILES raises ValueError")

    # ── MACCS fingerprint ──────────────────────────────────────────────────
    fp_maccs = mol.maccs()
    check(isinstance(fp_maccs, bytes), "maccs() returns bytes")
    check(len(fp_maccs) == 21, f"maccs() returns 21 bytes (got {len(fp_maccs)})")

    # benzene MACCS should have bit 162 (aromatic ring key) set
    byte_idx, bit_off = 162 // 8, 162 % 8
    check(bool(fp_maccs[byte_idx] & (1 << bit_off)), "benzene has MACCS bit 162 (aromatic)")

    # ── Morgan fingerprint ─────────────────────────────────────────────────
    fp_morgan = mol.morgan(radius=2, nbits=2048)
    check(isinstance(fp_morgan, bytes), "morgan() returns bytes")
    check(len(fp_morgan) == 256, f"morgan(nbits=2048) returns 256 bytes (got {len(fp_morgan)})")
    check(any(b != 0 for b in fp_morgan), "morgan fingerprint is non-zero for benzene")

    fp_morgan_512 = mol.morgan(radius=2, nbits=512)
    check(len(fp_morgan_512) == 64, f"morgan(nbits=512) returns 64 bytes (got {len(fp_morgan_512)})")

    # ── tanimoto ──────────────────────────────────────────────────────────
    sim_self = molprint.tanimoto(fp_maccs, fp_maccs)
    check(abs(sim_self - 1.0) < 1e-9, f"tanimoto(x, x) == 1.0 (got {sim_self})")

    mol2 = molprint.Mol.from_smiles("CCO")
    fp2 = mol2.maccs()
    sim = molprint.tanimoto(fp_maccs, fp2)
    check(0.0 <= sim < 1.0, f"tanimoto(benzene, ethanol) in (0, 1) (got {sim:.4f})")

    try:
        molprint.tanimoto(fp_maccs, fp_morgan)
        check(False, "mismatched lengths should raise ValueError")
    except ValueError:
        check(True, "mismatched fingerprint lengths raise ValueError")

    # ── batch_maccs ────────────────────────────────────────────────────────
    smiles = ["c1ccccc1", "CCO", "not_valid", "CC(=O)O"]
    results = molprint.batch_maccs(smiles)
    check(len(results) == 4, f"batch_maccs returns list of length 4 (got {len(results)})")
    check(results[0] is not None and len(results[0]) == 21, "batch_maccs[0] is 21 bytes")
    check(results[2] is None, "batch_maccs skips invalid SMILES (returns None)")

    # ── batch_morgan ───────────────────────────────────────────────────────
    morgan_results = molprint.batch_morgan(["c1ccccc1", "CCO"], radius=2, nbits=1024)
    check(len(morgan_results) == 2, "batch_morgan returns 2 results")
    check(all(r is not None and len(r) == 128 for r in morgan_results),
          "batch_morgan(nbits=1024) gives 128 bytes each")

    print()
    print("All tests passed.")


if __name__ == "__main__":
    main()
