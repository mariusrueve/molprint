# Fingerprints

molprint implements two fingerprint types: Morgan (circular/ECFP) and MACCS-166 structural keys.

## `FingerprintBits`

All fingerprints are represented as `FingerprintBits`, a word-aligned bit vector:

```rust
pub struct FingerprintBits {
    words: Vec<u64>,
    nbits: usize,
}
```

Operations are implemented at the word level using Rust's `u64::count_ones()` (which compiles to a single `POPCNT` instruction on x86/ARM). Bitwise AND and OR iterate over word pairs, making Tanimoto computation efficient for any bit width.

The hex encoding used in FPS files stores bytes least-significant first, matching chemfp's convention.

## The `Fingerprinter` trait

```rust
pub trait Fingerprinter {
    fn fingerprint(&self, mol: &MolGraph) -> FingerprintBits;
    fn name(&self) -> &str;
    fn nbits(&self) -> usize;
}
```

Both `Morgan` and `Maccs166` implement this trait. The CLI uses `Box<dyn Fingerprinter>` to select the algorithm at runtime.

## Morgan / ECFP

Morgan fingerprints encode circular atom environments. For each atom, an environment of radius `r` includes everything within `r` bonds.

### Atom invariants (radius 0)

Each atom gets an initial hash from its local properties:

```rust
fn atom_invariant(&self, mol: &MolGraph, idx: NodeIndex) -> u32 {
    let atom = mol.atom(idx);
    let packed = ((atom.element.atomic_number() as u64) << 48)
        | ((mol.degree(idx) as u64) << 32)
        | ((atom.h_count as u64) << 16)
        | (((atom.charge as i64 + 128) as u64) << 8)
        | (atom.in_ring as u64);
    hash_u64(packed) as u32
}
```

The properties encoded are: atomic number, heavy-atom degree, hydrogen count, formal charge, ring membership. This matches the standard ECFP invariants.

### Iteration

At each radius step, each atom's hash is updated based on its neighbors:

```rust
fn iterate(&self, mol: &MolGraph, prev_hashes: &[u32], _radius: u32) -> Vec<u32> {
    for i in 0..n {
        let mut neighbors: Vec<(u8, u32)> = mol.edges(idx)
            .map(|e| (bond_order, prev_hashes[neighbor]))
            .collect();
        neighbors.sort_unstable(); // deterministic ordering
        new_hashes[i] = combine_hashes(prev_hashes[i] as u64, &neighbors);
    }
}
```

Sorting neighbors before hashing is critical for canonical fingerprints — without it, the result would depend on the order atoms were added to the graph.

### Hash function

The hash uses a FNV-1a-inspired 64-bit hash for atom invariants, and boost-style `hash_combine` for neighborhood accumulation:

```rust
fn hash_combine(h: u64, val: u64) -> u64 {
    h ^ val.wrapping_add(0x9e3779b97f4a7c15)
         .wrapping_add(h << 6)
         .wrapping_add(h >> 2)
}
```

**Important**: `AHasher::default()` uses a randomized seed and cannot be used here. All hashing must be deterministic across process launches.

### Bit folding

All hashes from radius 0 through `radius` are collected into a `HashSet<u32>`, then each is folded into the bit vector:

```rust
for h in all_hashes {
    fp.set(h as usize % self.nbits);
}
```

Using a `HashSet` means each unique environment contributes at most one bit position, even if it appears at multiple atoms. This is equivalent to RDKit's default behavior.

### Configurations

| Name | `Morgan::new(r, bits)` |
|---|---|
| ECFP2 | `Morgan::new(1, 2048)` |
| ECFP4 | `Morgan::new(2, 2048)` |
| ECFP6 | `Morgan::new(3, 2048)` |

## MACCS-166

MACCS-166 is a fixed set of 166 structural keys defined as patterns (originally as SMARTS queries in RDKit). Each key tests for the presence of a specific structural feature.

The keys are implemented manually in `maccs.rs` rather than via runtime SMARTS evaluation. This provides:

1. Exact control over counting semantics (`uniquify=True, min_count=1` matching RDKit)
2. Better performance (no SMARTS compilation overhead per fingerprint)
3. Direct bit-level RDKit compatibility

### Implementation structure

Each bit is computed by a dedicated test. The tests fall into several categories:

- **Element presence**: "does the molecule contain a sulfur atom?"
- **Ring tests**: "is there a 6-membered ring?", "is there an aromatic ring?"
- **Functional group patterns**: carbonyl, hydroxyl, amine, nitro, etc.
- **Chain patterns**: specific bond sequences, e.g., N-C=O (amide)
- **Count-based keys**: "more than 2 ring systems", etc.

A helper function `has_simple_path` supports keys that require checking for specific atom sequences:

```rust
fn has_simple_path(
    mol: &MolGraph,
    start: NodeIndex,
    length: usize,
    target_pred: &dyn Fn(NodeIndex) -> bool,
) -> bool
```

For molecules with ≤ 64 atoms (the vast majority), this uses a single `u64` bitmask for the visited set, avoiding heap allocation.

### Accuracy

MACCS-166 achieves 100% accuracy on the ChEMBL 10k validation set against RDKit. Getting there required careful attention to:

- **uniquify=True**: RDKit's default MACCS implementation uses `uniquify=True` in its SMARTS matching, meaning each atom in the molecule is only used once per match. Some keys set the bit if the count is ≥ some threshold.
- **Aromatic implicit H**: aromatic nitrogen (pyridine vs. pyrrole) has different H counts, and several MACCS keys test for `[nH]` specifically.
- **Bit 0**: always 0 (undefined in the original MACCS spec). The bitvector has 167 bits, with bit 0 unused.

### Fingerprinter implementation

```rust
impl Fingerprinter for Maccs166 {
    fn nbits(&self) -> usize { 167 }
    fn name(&self) -> &str { "MACCS166" }
    fn fingerprint(&self, mol: &MolGraph) -> FingerprintBits {
        let mut fp = FingerprintBits::new(167);
        // bit 1: isotope
        // bit 2: ...
        // ...
        fp
    }
}
```
