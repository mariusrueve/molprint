# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — 2026-03-21

### Added

#### molprint-core
- `MolGraph` molecular graph type (petgraph `UnGraph` wrapper)
- `Atom` struct with element, charge, isotope, H count, aromaticity, ring membership
- `Element` enum covering H–Bi plus `Unknown`
- `BondType` enum: `Single`, `Double`, `Triple`, `Aromatic`
- SMILES lexer + parser: organic subset, bracket atoms, ring closures, branches,
  stereo bond tokens (ignored), multi-component (`.`)
- Ring perception (SSSR) via BFS-based Horton algorithm; circuit rank formula
- Aromaticity perception module
- **SMARTS v1**: compile SMARTS patterns, `has_match`, `subgraph_match`
  - Atom primitives: `#n`, element symbols, `a`/`A`, `Hn`, `+n`/`-n`,
    `R`/`R0`/`Rn`, `Dn`, `*`, logical `!`, `&`, `;`, `,`
  - Bond primitives: `-`, `=`, `#`, `:`, `~`, unspecified (single/aromatic)
  - VF2-style recursive subgraph isomorphism

#### molprint-fp
- `FingerprintBits`: word-aligned bit vector with POPCNT, hex serialisation
- `Fingerprinter` trait for composable fingerprint algorithms
- `Morgan`: iterative environment hashing (configurable radius + bit count)
- `Maccs166`: 166-bit MACCS structural keys, validated bit-exact against RDKit

#### molprint-search
- `tanimoto`, `dice`, `cosine` similarity functions (u64 POPCNT loops)
- `threshold_search` and `top_k_search` via rayon parallel iteration

#### molprint-io
- `SmilesFileReader`: streaming iterator for SMILES files (`SMILES\tID`)
- `SdfReader`: streaming SDF V2000 parser with gzip (`.sdf.gz`) support
- `FpsReader` / `FpsWriter`: chemfp-compatible FPS format

#### molprint-cli
- `molprint fp`: compute fingerprints from SMILES/SDF → FPS
- `molprint search`: query SMILES + FPS database → ranked similarity hits

#### Python bindings (molprint-py)
- `Mol.from_smiles(smiles)` — parse SMILES, raises `ValueError` on error
- `mol.maccs()` → `bytes` (21 bytes, 166 bits)
- `mol.morgan(radius, nbits)` → `bytes`
- `tanimoto(fp_a, fp_b)` → `float`
- `batch_maccs(smiles_list)`, `batch_morgan(smiles_list, radius, nbits)`
- `pip install .` via maturin + `pyproject.toml`

#### Infrastructure
- Criterion benchmarks: SMILES parsing, Morgan/MACCS throughput, Tanimoto
- cargo-fuzz target for SMILES parser (`fuzz/fuzz_targets/fuzz_smiles.rs`)
- CI: test, clippy, fmt, doc, MSRV (1.75), security audit
- Scheduled CI: fuzz (60 s), RDKit MACCS validation
- Scripts: ChEMBL corpus download, corpus validation, RDKit benchmark comparison
