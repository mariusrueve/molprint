# Similarity & Screening

The `molprint-search` crate provides similarity metrics and parallel screening functions.

## Similarity metrics

All three metrics operate on `FingerprintBits` and use POPCNT for efficiency.

### Tanimoto (Jaccard)

The standard similarity metric in cheminformatics:

```
T(A, B) = |A ∩ B| / |A ∪ B| = popcount(A & B) / popcount(A | B)
```

```rust
pub fn tanimoto(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    let intersection = a.and(b).count_ones() as f64;
    let union = a.or(b).count_ones() as f64;
    if union == 0.0 { 0.0 } else { intersection / union }
}
```

Returns 0.0 for two all-zero fingerprints (rather than NaN). Returns 1.0 when the fingerprints are identical.

A common threshold for "similar" compounds is Tanimoto ≥ 0.7 for ECFP4.

### Dice

```
D(A, B) = 2 * |A ∩ B| / (|A| + |B|)
```

Dice is more lenient than Tanimoto for compounds of different sizes, because the denominator uses the sum of individual set sizes rather than the union.

### Cosine

```
C(A, B) = |A ∩ B| / sqrt(|A| * |B|)
```

Cosine similarity treats the fingerprint as a vector and measures the angle between them. It is less commonly used in cheminformatics than Tanimoto but useful for some applications.

### Performance

At 36 ns per Tanimoto comparison (2048-bit), the bottleneck is memory bandwidth, not computation. The `u64` word loop uses Rust's `u64::count_ones()` which compiles to a single `POPCNT` instruction. For future work, explicit SIMD (AVX2 or NEON) could batch multiple fingerprint comparisons per loop iteration.

## Parallel screening

### `threshold_search`

Returns all database fingerprints with Tanimoto ≥ threshold, sorted by descending similarity:

```rust
pub fn threshold_search(
    query: &FingerprintBits,
    database: &[FingerprintBits],
    threshold: f64,
) -> Vec<SearchHit>
```

Uses Rayon's `par_iter()` for parallel evaluation. Each thread processes a chunk of the database independently, then results are merged and sorted.

### `top_k_search`

Returns the k most similar fingerprints:

```rust
pub fn top_k_search(
    query: &FingerprintBits,
    database: &[FingerprintBits],
    k: usize,
) -> Vec<SearchHit>
```

Computes similarity for every database entry (in parallel), sorts, then truncates to k. For very large databases where k is much smaller than the total size, a heap-based algorithm would be more efficient — this is the simplest correct implementation.

### `SearchHit`

```rust
pub struct SearchHit {
    pub index: usize,      // index into the database slice
    pub similarity: f64,
}
```

The index into the original `database` slice lets the caller look up the molecule ID.

## Screening performance

On Apple M-series, screening a 100k compound database takes ~826 µs per query (121M fingerprint comparisons/second). This is primarily memory-bound. At 2048 bits = 256 bytes per fingerprint, a 100k database is 25 MB — fits in L3 cache on most modern CPUs, which explains why Rayon parallelism helps.

For databases that fit in RAM (up to ~10M compounds at 2048-bit), molprint keeps everything in a `Vec<FingerprintBits>` without indexing structures. For larger scales, pre-filtering using popcount bounds (the Tanimoto lower bound from bit count differences) would reduce the number of full comparisons needed.

## Usage example

```rust
use molprint_search::screen::{threshold_search, top_k_search};

// database is Vec<FingerprintBits>, loaded from an FPS file
// query_fp is FingerprintBits for the query molecule

// All compounds with Tanimoto >= 0.6
let hits = threshold_search(&query_fp, &database, 0.6);
for hit in &hits {
    println!("index={} sim={:.4}", hit.index, hit.similarity);
}

// Top 10 most similar
let top10 = top_k_search(&query_fp, &database, 10);
```
