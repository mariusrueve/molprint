use crate::mol::atom::Element;
use crate::mol::bond::BondType;
use crate::mol::graph::MolGraph;
use crate::ring::find_sssr;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

/// Perceive aromaticity using Hückel's rule (4n+2 π electrons).
///
/// For each ring in SSSR, count π electrons:
/// - C with a double bond to a ring neighbor: 1 electron
/// - C marked aromatic (from SMILES): 1 electron
/// - N with lone pair (in ring, not double bond to ring neighbor): 2 electrons
/// - N with double bond to ring neighbor: 1 electron
/// - O in ring: 2 electrons
/// - S in ring: 2 electrons
///
/// If total π electrons = 4n+2 for some n ≥ 0, the ring is aromatic.
pub fn perceive_aromaticity(mol: &mut MolGraph) {
    let rings = find_sssr(mol);
    for ring in &rings {
        let ring_set: std::collections::HashSet<NodeIndex> = ring.iter().copied().collect();
        let pi_electrons = count_pi_electrons(mol, ring, &ring_set);

        // Hückel: 4n+2 → aromatic
        let aromatic = matches!(pi_electrons, 2 | 6 | 10 | 14 | 18 | 22);

        if aromatic {
            let n = ring.len();
            for i in 0..n {
                let u = ring[i];
                let v = ring[(i + 1) % n];
                mol[u].aromatic = true;
                if let Some(edge) = mol.find_edge(u, v) {
                    mol[edge] = BondType::Aromatic;
                }
            }
        }
    }
}

fn count_pi_electrons(
    mol: &MolGraph,
    ring: &[NodeIndex],
    ring_set: &std::collections::HashSet<NodeIndex>,
) -> u32 {
    let mut pi = 0u32;

    for &atom_idx in ring {
        let atom = &mol[atom_idx];

        // Check if this atom has a double bond to a ring neighbor
        let has_double_to_ring = mol.edges(atom_idx).any(|e| {
            let nb = if e.source() == atom_idx {
                e.target()
            } else {
                e.source()
            };
            ring_set.contains(&nb) && *e.weight() == BondType::Double
        });

        // If atom is already marked aromatic by SMILES
        let already_aromatic = atom.aromatic;

        pi += match atom.element {
            Element::C => {
                if has_double_to_ring || already_aromatic {
                    1
                } else {
                    0
                }
            }
            Element::N => {
                if has_double_to_ring {
                    1
                } else {
                    2
                }
            }
            Element::O | Element::S => 2,
            _ => {
                if has_double_to_ring || already_aromatic {
                    1
                } else {
                    0
                }
            }
        };
    }

    pi
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn test_benzene_aromatic() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        perceive_aromaticity(&mut mol);
        for node in mol.node_indices() {
            assert!(
                mol[node].aromatic,
                "atom {} should be aromatic",
                node.index()
            );
        }
    }

    #[test]
    fn test_pyridine_aromatic() {
        let mut mol = parse_smiles("c1ccncc1").unwrap();
        perceive_aromaticity(&mut mol);
        for node in mol.node_indices() {
            assert!(mol[node].aromatic);
        }
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        let mut mol = parse_smiles("C1CCCCC1").unwrap();
        perceive_aromaticity(&mut mol);
        for node in mol.node_indices() {
            assert!(!mol[node].aromatic);
        }
    }

    #[test]
    fn test_furan_aromatic() {
        let mut mol = parse_smiles("c1ccoc1").unwrap();
        perceive_aromaticity(&mut mol);
        for node in mol.node_indices() {
            assert!(mol[node].aromatic);
        }
    }
}
