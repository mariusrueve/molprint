# Ring Perception

Ring perception identifies which atoms belong to rings, what size those rings are, and produces the Smallest Set of Smallest Rings (SSSR). This information is used directly by the MACCS fingerprint (many keys are ring-related) and by the Morgan fingerprint (atom invariants include `in_ring`).

## Why SSSR?

The SSSR is the minimal basis for all rings in a molecule. For a molecule with `e` edges, `n` nodes, and `c` connected components, the circuit rank (number of independent rings) is:

```
r = e − n + c
```

For benzene: `6 − 6 + 1 = 1`. For naphthalene: `11 − 10 + 1 = 2`. The SSSR contains exactly `r` rings.

Alternative approaches (DFS cycle enumeration, relevant cycles) enumerate more rings. For naphthalene, they would find 3 rings (the two 6-membered rings plus the 10-membered outer ring). The SSSR approach matches RDKit's behavior, which is important for MACCS accuracy.

## Algorithm

The implementation uses a BFS-based approach sometimes called the Horton algorithm:

**Step 1: Generate candidates**

For each edge `(u, v)`:
1. Temporarily remove the edge
2. Find the shortest path from `u` to `v` using BFS (ignoring the removed edge)
3. If a path exists, the path + the removed edge forms a ring candidate

This produces at most one candidate per edge. Candidates are sorted by ring size (smallest first).

**Step 2: Select a basis**

Greedily select candidates to form a linearly independent set in GF(2) (the field with two elements). Each ring is represented as a bit vector over edges (1 if the ring contains that edge, 0 otherwise). A candidate is added to the SSSR if its edge vector is linearly independent from all previously selected rings, determined by Gaussian elimination over GF(2).

Selection stops when `r` rings have been chosen.

## Implementation details

```rust
pub fn find_sssr(mol: &MolGraph) -> Vec<Vec<NodeIndex>> {
    let target_rings = e.saturating_sub(n) + components;
    // ... generate candidates via BFS ...
    // ... select basis via GF(2) Gaussian elimination ...
}
```

The GF(2) basis is maintained in row-echelon form. `reduce_to_basis` XORs the candidate vector against basis vectors whose pivot column is set in the candidate. If the result is non-zero, the candidate is independent and is added.

## Ring info assignment

After `find_sssr`, the `assign_ring_info` function marks atoms:

```rust
pub fn assign_ring_info(mol: &mut MolGraph) {
    for ring in &rings {
        let ring_size = ring.len() as u8;
        for &atom_idx in ring {
            mol[atom_idx].in_ring = true;
            mol[atom_idx].ring_sizes.push(ring_size);
        }
    }
}
```

Junction atoms (atoms shared between two fused rings, like those in naphthalene) get multiple entries in `ring_sizes`. This allows MACCS keys to correctly identify fused ring systems.

## Spiro atom handling

Spiro compounds (two rings sharing exactly one atom) were an early bug source. A naive ring detection approach that excludes atoms already in one ring would fail to detect the second ring in a spiro compound. The BFS candidate generation correctly handles this because it operates on edges, not atoms — the spiro center is visited in candidates for both rings independently.

The test `test_spiro_compound` covers `C1CCC2(CC1)CCCC2` (spiro[5.5]undecane) and asserts `find_sssr` returns exactly 2 rings.

## Connected components

The circuit rank formula requires knowing the number of connected components. Components are counted by a simple BFS that starts from each unvisited node:

```rust
fn count_components(mol: &MolGraph) -> usize {
    let mut visited = vec![false; mol.node_count()];
    let mut count = 0;
    for start in mol.node_indices() {
        if !visited[start.index()] {
            count += 1;
            // BFS from start...
        }
    }
    count
}
```

For disconnected SMILES like `CC.OO`, components = 2, and the circuit rank is 0 (no rings), which is correct.
