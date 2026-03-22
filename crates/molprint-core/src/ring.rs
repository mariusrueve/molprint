use crate::mol::graph::MolGraph;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::{HashMap, VecDeque};

/// Find the Smallest Set of Smallest Rings (SSSR) using a BFS-based approach.
///
/// For each edge (u, v): temporarily remove it, find shortest path u→v via BFS.
/// If a path exists, that path + the removed edge forms a ring candidate.
/// Filter candidates to the minimum independent set (circuit rank = edges - nodes + components).
pub fn find_sssr(mol: &MolGraph) -> Vec<Vec<NodeIndex>> {
    if mol.edge_count() == 0 || mol.node_count() == 0 {
        return vec![];
    }

    let n = mol.node_count();
    let e = mol.edge_count();

    // Count connected components via BFS
    let components = count_components(mol);
    let target_rings = e.saturating_sub(n) + components;

    if target_rings == 0 {
        return vec![];
    }

    // Generate ring candidates: for each edge, BFS without that edge
    let mut candidates: Vec<Vec<NodeIndex>> = Vec::new();
    for edge in mol.edge_indices() {
        let (u, v) = mol.edge_endpoints(edge).unwrap();
        if let Some(path) = bfs_path_without_edge(mol, u, v, edge) {
            let mut ring = path;
            ring.push(u); // close the ring
            candidates.push(ring);
        }
    }

    // Sort by ring size (smallest first)
    candidates.sort_by_key(|r| r.len());

    // Greedily select linearly independent rings
    let mut sssr: Vec<Vec<NodeIndex>> = Vec::new();
    let mut basis: Vec<Vec<bool>> = Vec::new(); // edge-membership vectors

    let edge_idx: HashMap<(usize, usize), usize> = mol
        .edge_indices()
        .enumerate()
        .map(|(i, eid)| {
            let (u, v) = mol.edge_endpoints(eid).unwrap();
            let (a, b) = (u.index().min(v.index()), u.index().max(v.index()));
            ((a, b), i)
        })
        .collect();

    for ring in &candidates {
        if sssr.len() >= target_rings {
            break;
        }
        let vec = ring_to_edge_vector(ring, &edge_idx, e);
        if let Some(reduced) = reduce_to_basis(&vec, &basis) {
            basis.push(reduced);
            sssr.push(ring.clone());
        }
    }

    sssr
}

/// BFS from `start` to `end`, ignoring the edge `skip_edge`.
/// Returns path from `start` to `end` (not including `start` at the beginning,
/// but including `end`), or None if no path exists.
fn bfs_path_without_edge(
    mol: &MolGraph,
    start: NodeIndex,
    end: NodeIndex,
    skip_edge: petgraph::graph::EdgeIndex,
) -> Option<Vec<NodeIndex>> {
    let mut visited: HashMap<NodeIndex, NodeIndex> = HashMap::new();
    let mut queue = VecDeque::new();
    queue.push_back(start);
    visited.insert(start, start);

    while let Some(current) = queue.pop_front() {
        for edge_ref in mol.edges(current) {
            if edge_ref.id() == skip_edge {
                continue;
            }
            let neighbor = if edge_ref.source() == current {
                edge_ref.target()
            } else {
                edge_ref.source()
            };

            if visited.contains_key(&neighbor) {
                continue;
            }
            visited.insert(neighbor, current);

            if neighbor == end {
                // Reconstruct path from end to start
                let mut path = vec![end];
                let mut node = end;
                while node != start {
                    node = visited[&node];
                    if node != start {
                        path.push(node);
                    }
                }
                path.reverse();
                return Some(path);
            }
            queue.push_back(neighbor);
        }
    }
    None
}

fn count_components(mol: &MolGraph) -> usize {
    let mut visited = vec![false; mol.node_count()];
    let mut count = 0;
    for start in mol.node_indices() {
        if !visited[start.index()] {
            count += 1;
            let mut queue = VecDeque::new();
            queue.push_back(start);
            visited[start.index()] = true;
            while let Some(node) = queue.pop_front() {
                for nb in mol.neighbors(node) {
                    if !visited[nb.index()] {
                        visited[nb.index()] = true;
                        queue.push_back(nb);
                    }
                }
            }
        }
    }
    count
}

/// Convert a ring (as a list of atoms) to a bit vector over edges.
fn ring_to_edge_vector(
    ring: &[NodeIndex],
    edge_idx: &HashMap<(usize, usize), usize>,
    num_edges: usize,
) -> Vec<bool> {
    let mut vec = vec![false; num_edges];
    let n = ring.len();
    for i in 0..n {
        let u = ring[i].index();
        let v = ring[(i + 1) % n].index();
        let key = (u.min(v), u.max(v));
        if let Some(&idx) = edge_idx.get(&key) {
            vec[idx] = true;
        }
    }
    vec
}

/// Reduce `vec` against the current GF(2) row-echelon basis.
///
/// Returns the reduced vector (non-zero) if `vec` is linearly independent from
/// the basis, or `None` if it reduces to the zero vector (dependent).
/// The basis must be in row-echelon form (each vector's first set bit is
/// unique); pushing the returned vector maintains that invariant.
fn reduce_to_basis(vec: &[bool], basis: &[Vec<bool>]) -> Option<Vec<bool>> {
    let mut v = vec.to_vec();
    for b in basis {
        let pivot = b.iter().position(|&x| x);
        if let Some(p) = pivot {
            if v[p] {
                for i in 0..v.len() {
                    v[i] ^= b[i];
                }
            }
        }
    }
    if v.iter().any(|&x| x) {
        Some(v)
    } else {
        None
    }
}

/// Mark atoms and ring sizes on the molecular graph.
pub fn assign_ring_info(mol: &mut MolGraph) {
    let rings = find_sssr(mol);
    for node in mol.node_indices() {
        mol[node].in_ring = false;
        mol[node].ring_sizes.clear();
    }
    for ring in &rings {
        let ring_size = ring.len() as u8;
        for &atom_idx in ring {
            mol[atom_idx].in_ring = true;
            // Allow duplicates so junction atoms record each ring they belong to,
            // enabling detection of same-size fused rings (e.g., naphthalene 6+6).
            mol[atom_idx].ring_sizes.push(ring_size);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn ring_debug_fn_molecules() {
        let smiles = [
            ("COC(=O)CC1NN=C2N(CCN2c2cccc(Cl)c2)C1=O", "mol1"),
            ("CCCC1(CC(=O)O)CCCc2c1[nH]c1c(F)ccc(C#N)c21", "mol2"),
        ];
        for (smi, name) in &smiles {
            let mol = parse_smiles(smi).unwrap();
            let n_atoms = mol.node_count();
            let n_in_ring = mol.node_indices().filter(|&n| mol[n].in_ring).count();
            let mut all_ring_sizes: Vec<u8> = mol
                .node_indices()
                .flat_map(|n| mol[n].ring_sizes.iter().copied())
                .collect::<std::collections::HashSet<_>>()
                .into_iter()
                .collect();
            all_ring_sizes.sort();
            let sssr = find_sssr(&mol);
            eprintln!(
                "{name}: {n_in_ring}/{n_atoms} in_ring, SSSR={} rings, sizes={all_ring_sizes:?}",
                sssr.len()
            );
            for (i, ring) in sssr.iter().enumerate() {
                let indices: Vec<_> = ring.iter().map(|n| n.index()).collect();
                eprintln!("  Ring {i} (size {}): {indices:?}", ring.len());
            }
            // Print all atoms and their neighbor indices
            eprintln!("  Atom list:");
            for n in mol.node_indices() {
                let elem = format!("{:?}", mol[n].element);
                let neighbors: Vec<usize> = mol.neighbors(n).map(|nb| nb.index()).collect();
                eprintln!("    atom {}: {} neighbors={:?}", n.index(), elem, neighbors);
            }
        }
    }

    #[test]
    fn test_benzene_ring() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        assign_ring_info(&mut mol);
        assert_eq!(find_sssr(&mol).len(), 1);
        for node in mol.node_indices() {
            assert!(mol[node].in_ring);
            assert!(mol[node].ring_sizes.contains(&6));
        }
    }

    #[test]
    fn test_naphthalene() {
        let mut mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let rings = find_sssr(&mol);
        assign_ring_info(&mut mol);
        assert_eq!(rings.len(), 2);
    }

    #[test]
    fn test_no_rings() {
        let mol = parse_smiles("CCCC").unwrap();
        assert_eq!(find_sssr(&mol).len(), 0);
    }

    #[test]
    fn test_cyclopropane() {
        let mol = parse_smiles("C1CC1").unwrap();
        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 3);
    }

    #[test]
    fn test_spiro_compound() {
        let mol = parse_smiles("C1CCC2(CC1)CCCC2").unwrap();
        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 2);
    }
}
