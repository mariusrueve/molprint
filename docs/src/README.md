# molprint

**molprint** is a high-performance molecular fingerprint computation and similarity search library written in Rust.

The goal is to be fast enough for large-scale virtual screening — targeting 5–10× faster than RDKit for bulk workflows — while remaining accurate enough to match RDKit bit-for-bit on standard benchmarks.

## What it does

- Parses SMILES strings into molecular graphs
- Computes **MACCS-166** structural key fingerprints (100% RDKit-accurate on ChEMBL 10k)
- Computes **Morgan/ECFP** circular fingerprints at configurable radius and bit width (512–4096)
- Calculates **Tanimoto, Dice, and Cosine** similarity using POPCNT on `u64` word arrays
- Runs **parallel threshold and top-k screening** via Rayon
- Reads and writes **FPS, SDF, and SMILES** file formats (FPS is chemfp-compatible)

## Benchmarks

Measured on Apple M-series, Rust 1.94, `--release`.

| Operation | Performance |
|---|---|
| Tanimoto (2048-bit) | 36 ns |
| Morgan ECFP4 batch | ~700k mol/s |
| MACCS-166 batch | ~535k mol/s |
| Screening 100k compounds | 826 µs / query |

## Quick example

```rust
use molprint_core::smiles::parse_smiles;
use molprint_fp::{morgan::Morgan, maccs::Maccs166, traits::Fingerprinter};
use molprint_search::metrics::tanimoto;

let mol_a = parse_smiles("c1ccccc1").unwrap();   // benzene
let mol_b = parse_smiles("c1ccncc1").unwrap();   // pyridine

let fp = Morgan::new(2, 2048);
let sim = tanimoto(&fp.fingerprint(&mol_a), &fp.fingerprint(&mol_b));
println!("{:.3}", sim); // 0.600
```

## Workspace layout

| Crate | Role |
|---|---|
| `molprint-core` | SMILES parser, molecular graph, ring perception, SMARTS |
| `molprint-fp` | Morgan and MACCS fingerprint algorithms |
| `molprint-search` | Similarity metrics and parallel screening |
| `molprint-io` | FPS, SDF, SMILES file I/O |
| `molprint-cli` | End-user CLI tool |

Crates form a strict dependency chain: `molprint-core` → `molprint-fp` → `molprint-search` + `molprint-io` → `molprint-cli`. No cycles, no dev-dependency shortcuts across this chain.

## License

MIT
