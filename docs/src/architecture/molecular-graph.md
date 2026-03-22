# Molecular Graph

The molecular graph is the central data structure in molprint. Every algorithm — ring perception, fingerprinting, SMARTS matching — operates on it.

## Type definition

```rust
// mol/graph.rs
pub type MolGraph = UnGraph<Atom, BondType>;
```

`MolGraph` is a type alias for petgraph's `UnGraph<Atom, BondType>`, an undirected graph where each node carries an `Atom` and each edge carries a `BondType`. Hydrogen atoms are not stored as nodes in the graph (they are implicit); instead, `Atom::h_count` records the total hydrogen count.

## The `Atom` struct

```rust
pub struct Atom {
    pub element: Element,    // atomic number-backed enum
    pub charge: i8,          // formal charge
    pub isotope: Option<u16>,
    pub aromatic: bool,
    pub h_count: u8,         // total H count (implicit + explicit)
    pub explicit_h: Option<u8>, // Some if bracket atom with H spec
    pub map_num: u16,        // atom-map number (for reaction SMILES)
    pub in_ring: bool,       // set by ring perception
    pub ring_sizes: Vec<u8>, // one entry per SSSR ring containing this atom
}
```

`h_count` is computed after the full graph is constructed by `assign_implicit_hydrogens()`. For bracket atoms (e.g., `[NH3]`), the H count is explicit; for organic subset atoms (e.g., `N` in `CC(N)C`), it is inferred from the atom's default valence minus its current bond order sum.

## Bond types

```rust
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
}
```

Each `BondType` has a `valence_contribution()` method that returns its integer bond order (Single=1, Double=2, Triple=3, Aromatic=1). The aromatic contribution of 1 is intentional: the additional electron from the pi system is accounted for separately in the implicit H calculation.

## The `MolGraphExt` trait

petgraph's `UnGraph` is a generic graph; it doesn't know about chemistry. The `MolGraphExt` trait adds chemistry-aware methods:

```rust
pub trait MolGraphExt {
    fn num_atoms(&self) -> usize;
    fn num_bonds(&self) -> usize;
    fn atom(&self, idx: NodeIndex) -> &Atom;
    fn atom_mut(&mut self, idx: NodeIndex) -> &mut Atom;
    fn bond_between(&self, a: NodeIndex, b: NodeIndex) -> Option<BondType>;
    fn heavy_neighbors(&self, idx: NodeIndex) -> Vec<NodeIndex>;
    fn degree(&self, idx: NodeIndex) -> usize;
    fn compute_implicit_h(&self, idx: NodeIndex) -> u8;
    fn assign_implicit_hydrogens(&mut self);
}
```

`compute_implicit_h` is the most subtle method. It handles three cases:

1. **Bracket atom with H spec** (e.g., `[NH3]`): return the explicit count directly.
2. **Non-organic-subset element** (e.g., `[Fe]`): return 0. These atoms get no implicit Hs.
3. **Organic subset atom**: find the smallest standard valence ≥ current bond order sum; return `valence − bond_sum − |charge|`. For aromatic atoms, an extra +1 is added to the effective bond sum to account for the pi-electron contribution.

## Implicit hydrogen calculation for aromatic atoms

This was the single most difficult piece to get right. The logic is:

```
effective_bond_sum = bond_sum + 1   (add 1 for the pi electron)
h_count = valence - effective_bond_sum
```

But this +1 must only be applied when the atom's valence is not already fully saturated by its sigma bonds. For example, trimethylamine-oxide's nitrogen in an aromatic context is a different case than a pyridine nitrogen. The code checks whether `valence >= bond_sum` before deciding to apply the adjustment. If the valence is already fully consumed by heavy-atom bonds, no pi adjustment is added, and 0 implicit Hs are assigned.

## Node and edge indexing

petgraph assigns `NodeIndex` and `EdgeIndex` values sequentially as nodes and edges are added. These indices are stable within a single `MolGraph` lifetime (no removals happen in molprint). The SMILES parser, ring perception, and all fingerprint algorithms use `NodeIndex::new(i)` to iterate over atoms.
