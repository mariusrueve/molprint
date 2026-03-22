# molprint

Molecular fingerprint computation and similarity search in Rust — fast enough for large-scale screening, accurate enough to match RDKit bit-for-bit.

- **MACCS-166** structural keys (100% RDKit-accurate on ChEMBL 10k)
- **Morgan/ECFP** fingerprints (configurable radius, 512–4096 bit)
- **Tanimoto, Dice, Cosine** similarity via POPCNT
- **Parallel threshold and top-k screening** via Rayon
- **FPS, SDF, SMILES** file I/O (chemfp-compatible)

→ **[Documentation](https://mariusrueve.github.io/molprint/)**

## Install

**From a release binary** — download from [Releases](https://github.com/mariusrueve/molprint/releases), unpack, and place `molprint` on your `PATH`.

> **macOS note:** binaries are ad-hoc signed but not notarized. After downloading, clear the quarantine flag once:
> ```bash
> xattr -dr com.apple.quarantine molprint
> ```

**From source:**

```bash
cargo install --git https://github.com/mariusrueve/molprint molprint-cli
```

## CLI

```bash
# Compute fingerprints
molprint fp --input molecules.smi --output molecules.fps
molprint fp --input molecules.smi --fp-type maccs --output maccs.fps
molprint fp --input molecules.sdf --radius 3 --nbits 4096 --output ecfp6.fps

# Similarity search
molprint search --query "c1ccccc1" --db molecules.fps --threshold 0.4
molprint search --query "c1ccccc1" --db molecules.fps --top-k 20
```

Output: tab-separated `id\tsimilarity`, sorted by descending score.

## Library

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

## Benchmarks

Measured on Apple M-series, Rust 1.94, `--release`. Run `cargo bench` to reproduce.

| | |
|---|---|
| Tanimoto (2048-bit) | 36 ns |
| Morgan ECFP4 | ~700k mol/s |
| MACCS-166 | ~535k mol/s |
| Screening 100k compounds | 826 µs / query |

## Workspace

| Crate | Role |
|---|---|
| `molprint-core` | SMILES parser, molecular graph, ring perception, SMARTS |
| `molprint-fp` | Morgan and MACCS fingerprint algorithms |
| `molprint-search` | Similarity metrics and parallel screening |
| `molprint-io` | FPS, SDF, SMILES file I/O |
| `molprint-cli` | CLI tool |

## License

MIT
