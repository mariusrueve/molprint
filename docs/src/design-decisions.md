# Design Decisions

This page explains the key architectural choices made in molprint and the reasoning behind them.

## petgraph vs. a custom graph

The molecular graph uses petgraph's `UnGraph` rather than a custom adjacency list or matrix.

**Why petgraph**: Writing a correct, fast, well-tested graph data structure is a large project in itself. petgraph is mature, widely used in the Rust ecosystem, and provides exactly the operations we need: node/edge iteration, neighbor traversal, edge lookup by endpoint pair, and `NodeIndex` stability. The BFS in ring perception, the VF2 matching in SMARTS, and the neighborhood iteration in Morgan all map cleanly onto petgraph's API.

**The cost**: petgraph's `UnGraph` stores adjacency as a list of (source, target, edge) triples, which is not the most cache-friendly layout for large-scale operations. But for molecules (typically < 100 atoms), the in-memory layout is small enough to fit in cache regardless, and the performance difference vs. a custom flat structure is negligible.

**The alternative**: A flat array of `(u8, u8, BondType)` triples indexed by a CSR (compressed sparse row) adjacency list would be more cache-friendly for bulk operations. This would matter most for the SMARTS matcher traversing large molecules, but for typical drug-like molecules the difference is not measurable.

## Manual MACCS-166 vs. runtime SMARTS

MACCS-166 keys are defined in RDKit as SMARTS patterns. The obvious implementation is to compile those 166 SMARTS patterns at startup and evaluate each one per molecule.

molprint instead implements the 166 keys as hand-written Rust code.

**Why manual**: Three reasons.

1. **Exact RDKit semantics**: RDKit's MACCS implementation uses specific counting logic (`uniquify=True`, `min_count` thresholds per key) that is not the same as simply checking `has_match`. The counting semantics are easier to express precisely in code than to replicate through SMARTS configuration.

2. **Performance**: Compiling 166 SMARTS patterns and running VF2 matching for each molecule adds significant constant overhead per molecule. Manual code eliminates the pattern compilation step and allows simpler, more direct logic for common tests (element presence, ring sizes) without the general-purpose SMARTS machinery.

3. **Accuracy as first-class concern**: Because the keys are hand-written, each one can be individually tested against RDKit output. The validation test suite (`validate_against_rdkit.rs`) checks all 166 bits on a ChEMBL 10k subset. Catching a semantic error in a hand-written key is straightforward; tracking down why a SMARTS-based implementation differs from RDKit's behavior would be harder.

## `Vec<u64>` bitvector

Fingerprints are stored as `Vec<u64>` with a fixed bit count. The alternatives were:

- **`[u64; N]` const generic array**: Would require const generic parameters through the entire codebase and cannot represent runtime-configurable bit widths (Morgan supports 512–4096 bits).
- **`bitvec` crate**: A well-designed library, but adds a dependency and has more API surface than needed. The custom `FingerprintBits` has exactly the operations needed (set, get, count_ones, and/xor, hex encode/decode) with no extra overhead.
- **`Vec<u8>` or `Vec<u32>`**: `u64` matches the natural word size of modern CPUs, and `u64::count_ones()` compiles to a single `POPCNT` instruction. Using `u8` or `u32` would require more iterations per Tanimoto calculation.

## Rayon for parallelism

The screening functions use Rayon's `par_iter()` rather than manual thread management or async.

**Why Rayon**: `threshold_search` over a database of fingerprints is embarrassingly parallel — each comparison is independent. Rayon's work-stealing thread pool handles load balancing automatically, and the API change from `iter()` to `par_iter()` is minimal. Manual `std::thread` management would require partitioning the database into chunks, handling the results from each thread, and managing thread lifetimes — all complexity that Rayon handles for free.

**Why not async**: Similarity search is CPU-bound, not I/O-bound. Async would add complexity without benefiting throughput.

## `thiserror` for error types

All public error types use `thiserror::Error` derive macros.

**Why thiserror**: It generates `std::error::Error` and `Display` implementations from the error type definition, which would otherwise be boilerplate. The derive-based approach enforces consistent error message formatting and makes it easy to chain errors from lower-level crates (e.g., `LexerError` wraps into `ParseError` with `#[from]`).

**Alternative**: `anyhow` is common in applications but is designed for error propagation rather than typed error handling. Library code should use typed errors so callers can match on specific variants; `anyhow` boxes the error, losing the type.

## Deterministic fingerprints

Morgan's hash functions use fixed constants (`0x9e3779b97f4a7c15` from the golden ratio, and a fixed FNV-1a seed `0xcbf29ce484222325`). This is explicitly documented in the code with a comment: `AHasher::default() uses a randomized seed and must NOT be used here`.

**Why this matters**: If the hash function is randomized, the same SMILES produces a different fingerprint in each process invocation. This makes it impossible to compare fingerprints computed in different runs (e.g., a database built yesterday vs. a query today). Deterministic fingerprints are a correctness requirement, not just a performance preference.

## Implicit-only hydrogen representation

Hydrogen atoms are not stored as nodes in `MolGraph`. They are counted in `Atom::h_count` but not represented as graph nodes.

**Why**: Most cheminformatics algorithms (ring perception, fingerprinting, SMARTS matching) operate on the heavy-atom graph. Storing explicit H nodes would double the graph size for typical organic molecules, slow down neighbor iteration, and require special-casing Hs in every algorithm.

**The exception**: Bracket atoms like `[2H]` (deuterium) can appear in SMILES and do get added as nodes, because they carry isotope information that may be structurally meaningful. The `assign_implicit_hydrogens` function counts these explicit H neighbors and includes them in `h_count`.
