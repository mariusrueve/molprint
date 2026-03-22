# Architecture

This section describes the internal design of molprint — how each component works and why it was built the way it was.

## Crate map

```
molprint-core
    mol/atom.rs      — Element enum, Atom struct
    mol/bond.rs      — BondType enum
    mol/graph.rs     — MolGraph (petgraph UnGraph wrapper) + MolGraphExt trait
    smiles/lexer.rs  — tokenizer
    smiles/parser.rs — token stream → MolGraph
    ring.rs          — SSSR via BFS-based Horton algorithm
    arom.rs          — aromaticity perception
    smarts/          — SMARTS query language, VF2-style matching
        ast.rs
        lexer.rs
        matcher.rs

molprint-fp
    bitvec.rs        — FingerprintBits (Vec<u64> bit vector)
    traits.rs        — Fingerprinter trait
    morgan.rs        — Morgan/ECFP iterative hashing
    maccs.rs         — MACCS-166 structural keys

molprint-search
    metrics.rs       — Tanimoto, Dice, Cosine
    screen.rs        — threshold_search, top_k_search (Rayon)

molprint-io
    smiles_file.rs   — streaming SMILES line reader
    sdf.rs           — streaming SDF V2000 parser (plain + gzip)
    fps.rs           — chemfp FPS read/write

molprint-cli
    main.rs          — clap CLI: fp + search subcommands
```

## Data flow

For fingerprint computation:

```
SMILES/SDF file
    → molprint-io (SmilesFileReader / SdfReader)
    → parse_smiles → MolGraph
    → Fingerprinter::fingerprint → FingerprintBits
    → FpsWriter → .fps file
```

For similarity search:

```
.fps file → FpsReader → Vec<FingerprintBits>
query SMILES → parse_smiles → MolGraph → FingerprintBits
threshold_search / top_k_search (Rayon) → Vec<SearchHit>
```

## Design principles

- **Zero unsafe code** in library crates — relying on petgraph and the standard library for all unsafe operations.
- **No unwrap in library code** — all errors propagate through `Result` using `thiserror`.
- **Deterministic fingerprints** — hash functions use fixed seeds; no random state that varies between process launches.
- **Separation of concerns** — the molecular graph knows nothing about fingerprints; fingerprints know nothing about I/O.
