# Overview

molprint is structured as a Cargo workspace with five crates arranged in a strict dependency chain.

## Dependency order

```
molprint-core
    ↓
molprint-fp
    ↓          ↓
molprint-search  molprint-io
         ↓       ↓
         molprint-cli
```

This ordering enforces clean separation: the core graph library has no knowledge of fingerprints, fingerprint algorithms have no knowledge of file formats, and the CLI is the only place where everything comes together.

## molprint-core

The foundation of the library. Responsible for:

- Representing molecules as graphs (`MolGraph` = petgraph `UnGraph<Atom, BondType>`)
- Parsing SMILES strings into those graphs
- Computing ring membership (SSSR)
- Perceiving aromaticity
- Matching SMARTS patterns against molecule graphs

Nothing in `molprint-core` is fingerprint-specific.

## molprint-fp

Implements fingerprint algorithms on top of `MolGraph`. Contains:

- `FingerprintBits` — a word-aligned bit vector backed by `Vec<u64>`
- `Fingerprinter` trait — the common interface all algorithms implement
- `Morgan` — iterative neighborhood hashing (Morgan/ECFP)
- `Maccs166` — 166 manually implemented structural key tests

## molprint-search

Similarity metrics and parallel screening. Depends only on `molprint-fp` (for `FingerprintBits`). Contains:

- `tanimoto`, `dice`, `cosine` — bit-vector similarity metrics using POPCNT
- `threshold_search`, `top_k_search` — Rayon-parallelized screening functions

## molprint-io

File format support. Depends on `molprint-core` (for `MolGraph` and `parse_smiles`) and `molprint-fp` (for `FingerprintBits`). Contains:

- `SmilesFileReader` — streaming iterator over `(id, MolGraph)` pairs from SMILES files
- `SdfReader` — streaming iterator over SDF records (plain and gzip)
- `FpsWriter` / `FpsReader` — chemfp-compatible FPS format

## molprint-cli

The end-user binary. Ties everything together using `clap` for argument parsing. Has two subcommands:

- `fp` — read molecules, compute fingerprints, write FPS
- `search` — read an FPS database, compute a query fingerprint, run parallel search

## Testing strategy

- Unit tests live alongside the code they test, in `#[cfg(test)]` modules
- Integration tests live in `tests/` at the workspace root (`cli_integration.rs`)
- A fuzz target (`fuzz/fuzz_targets/fuzz_smiles.rs`) exercises the SMILES parser with arbitrary input
- Fingerprint accuracy is validated by `crates/molprint-fp/tests/validate_against_rdkit.rs`, which compares MACCS and Morgan output against RDKit on a ChEMBL subset
