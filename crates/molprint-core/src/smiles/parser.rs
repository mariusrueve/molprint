use super::lexer::{LexerError, Token};
use crate::mol::atom::Atom;
use crate::mol::bond::BondType;
use crate::mol::graph::{MolGraph, MolGraphExt};
use crate::ring::assign_ring_info;
use petgraph::graph::NodeIndex;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("lexer error: {0}")]
    Lexer(#[from] LexerError),
    #[error("unmatched ring closure {0}")]
    UnmatchedRing(u8),
    #[error("unmatched branch parenthesis")]
    UnmatchedBranch,
    #[error("empty SMILES string")]
    Empty,
    #[error("bond without target atom")]
    DanglingBond,
}

/// Parse a SMILES string into a molecular graph.
pub fn parse(smiles: &str) -> Result<MolGraph, ParseError> {
    let tokens = super::lexer::tokenize(smiles)?;
    if tokens.is_empty() {
        return Err(ParseError::Empty);
    }

    let mut graph = MolGraph::new_undirected();
    let mut current: Option<NodeIndex> = None;
    let mut branch_stack: Vec<NodeIndex> = Vec::new();
    let mut pending_bond: Option<BondType> = None;
    // ring_id → (atom_index, bond_type_if_specified)
    let mut ring_map: std::collections::HashMap<u8, (NodeIndex, Option<BondType>)> =
        std::collections::HashMap::new();

    for token in &tokens {
        match token {
            Token::Atom { element, aromatic } => {
                let mut atom = Atom::new(*element);
                atom.aromatic = *aromatic;
                let idx = graph.add_node(atom);

                if let Some(prev) = current {
                    let bond =
                        pending_bond
                            .take()
                            .unwrap_or(if *aromatic && graph[prev].aromatic {
                                BondType::Aromatic
                            } else {
                                BondType::Single
                            });
                    graph.add_edge(prev, idx, bond);
                }
                current = Some(idx);
            }

            Token::BracketAtom {
                isotope,
                element,
                aromatic,
                h_count,
                charge,
                map_num,
                chirality: _,
            } => {
                let mut atom = Atom::new(*element);
                atom.isotope = *isotope;
                atom.aromatic = *aromatic;
                // Bracket atoms without H specification have 0 implicit H (SMILES spec).
                // Use Some(0) so compute_implicit_h knows this is a bracket atom.
                atom.explicit_h = Some(h_count.unwrap_or(0));
                atom.charge = *charge;
                atom.map_num = *map_num;
                let idx = graph.add_node(atom);

                if let Some(prev) = current {
                    let bond =
                        pending_bond
                            .take()
                            .unwrap_or(if *aromatic && graph[prev].aromatic {
                                BondType::Aromatic
                            } else {
                                BondType::Single
                            });
                    graph.add_edge(prev, idx, bond);
                }
                current = Some(idx);
            }

            Token::Bond(bt) => {
                pending_bond = Some(*bt);
            }

            Token::BondStereo(_) => {
                // Stereo info stored separately; treat as single bond direction.
                // Don't override pending_bond.
            }

            Token::OpenBranch => {
                if let Some(cur) = current {
                    branch_stack.push(cur);
                }
            }

            Token::CloseBranch => {
                current = branch_stack.pop();
                if current.is_none() {
                    return Err(ParseError::UnmatchedBranch);
                }
                pending_bond = None;
            }

            Token::RingClosure(ring_id) => {
                let cur = current.ok_or(ParseError::DanglingBond)?;
                if let Some((other, ring_bond)) = ring_map.remove(ring_id) {
                    let bond = pending_bond.take().or(ring_bond).unwrap_or(
                        if graph[cur].aromatic && graph[other].aromatic {
                            BondType::Aromatic
                        } else {
                            BondType::Single
                        },
                    );
                    graph.add_edge(cur, other, bond);
                } else {
                    ring_map.insert(*ring_id, (cur, pending_bond.take()));
                }
            }

            Token::Dot => {
                current = None;
                pending_bond = None;
            }
        }
    }

    if let Some((&ring_id, _)) = ring_map.iter().next() {
        return Err(ParseError::UnmatchedRing(ring_id));
    }

    if !branch_stack.is_empty() {
        return Err(ParseError::UnmatchedBranch);
    }

    graph.assign_implicit_hydrogens();
    assign_ring_info(&mut graph);

    Ok(graph)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_methane() {
        let mol = parse("C").unwrap();
        assert_eq!(mol.node_count(), 1);
        assert_eq!(mol.edge_count(), 0);
        assert_eq!(mol[NodeIndex::new(0)].h_count, 4);
    }

    #[test]
    fn test_ethanol() {
        let mol = parse("CCO").unwrap();
        assert_eq!(mol.node_count(), 3);
        assert_eq!(mol.edge_count(), 2);
        assert_eq!(mol[NodeIndex::new(0)].h_count, 3);
        assert_eq!(mol[NodeIndex::new(1)].h_count, 2);
        assert_eq!(mol[NodeIndex::new(2)].h_count, 1);
    }

    #[test]
    fn test_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(mol.node_count(), 6);
        assert_eq!(mol.edge_count(), 6);
        for i in 0..6 {
            assert!(mol[NodeIndex::new(i)].aromatic);
            assert_eq!(mol[NodeIndex::new(i)].h_count, 1);
        }
    }

    #[test]
    fn test_acetic_acid() {
        let mol = parse("CC(=O)O").unwrap();
        assert_eq!(mol.node_count(), 4);
        assert_eq!(mol.edge_count(), 3);
        let c1 = NodeIndex::new(1);
        let o1 = NodeIndex::new(2);
        let edge = mol.find_edge(c1, o1).unwrap();
        assert_eq!(mol[edge], BondType::Double);
    }

    #[test]
    fn test_cyclopropane() {
        let mol = parse("C1CC1").unwrap();
        assert_eq!(mol.node_count(), 3);
        assert_eq!(mol.edge_count(), 3);
    }

    #[test]
    fn test_charged_ammonium() {
        let mol = parse("[NH4+]").unwrap();
        assert_eq!(mol.node_count(), 1);
        assert_eq!(mol[NodeIndex::new(0)].charge, 1);
        assert_eq!(mol[NodeIndex::new(0)].h_count, 4);
    }

    #[test]
    fn test_disconnected_fragments() {
        let mol = parse("CC.OO").unwrap();
        assert_eq!(mol.node_count(), 4);
        assert_eq!(mol.edge_count(), 2);
    }

    #[test]
    fn test_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        assert_eq!(mol.node_count(), 13);
        assert_eq!(mol.edge_count(), 13);
    }

    #[test]
    fn test_unmatched_ring_error() {
        assert!(parse("C1CC").is_err());
    }

    #[test]
    fn test_unmatched_branch_error() {
        assert!(parse("CC(O").is_err());
    }
}
