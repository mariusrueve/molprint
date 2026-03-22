use super::atom::Atom;
use super::bond::BondType;
use petgraph::graph::{NodeIndex, UnGraph};

/// The core molecular graph type.
pub type MolGraph = UnGraph<Atom, BondType>;

/// Extension trait adding chemistry-aware methods to MolGraph.
pub trait MolGraphExt {
    fn num_atoms(&self) -> usize;
    fn num_bonds(&self) -> usize;
    fn atom(&self, idx: NodeIndex) -> &Atom;
    fn atom_mut(&mut self, idx: NodeIndex) -> &mut Atom;
    fn bond_between(&self, a: NodeIndex, b: NodeIndex) -> Option<BondType>;
    fn heavy_neighbors(&self, idx: NodeIndex) -> Vec<NodeIndex>;
    fn degree(&self, idx: NodeIndex) -> usize;

    /// Compute implicit H count for an atom.
    /// - If explicit_h is set (bracket atom with H spec), return that.
    /// - For bracket atoms without H spec (explicit_h is None AND atom is bracket),
    ///   we return 0 — caller decides based on context.
    /// - For organic subset atoms, compute from default valences:
    ///   implicit_H = min_valence_geq_bond_sum - bond_sum
    fn compute_implicit_h(&self, idx: NodeIndex) -> u8;

    /// Set h_count on all atoms. Call after graph is fully constructed.
    fn assign_implicit_hydrogens(&mut self);
}

impl MolGraphExt for MolGraph {
    fn num_atoms(&self) -> usize {
        self.node_count()
    }

    fn num_bonds(&self) -> usize {
        self.edge_count()
    }

    fn atom(&self, idx: NodeIndex) -> &Atom {
        &self[idx]
    }

    fn atom_mut(&mut self, idx: NodeIndex) -> &mut Atom {
        &mut self[idx]
    }

    fn bond_between(&self, a: NodeIndex, b: NodeIndex) -> Option<BondType> {
        self.find_edge(a, b).map(|e| self[e])
    }

    fn heavy_neighbors(&self, idx: NodeIndex) -> Vec<NodeIndex> {
        self.neighbors(idx).collect()
    }

    fn degree(&self, idx: NodeIndex) -> usize {
        self.neighbors(idx).count()
    }

    fn compute_implicit_h(&self, idx: NodeIndex) -> u8 {
        let atom = &self[idx];

        // Bracket atoms with explicit H count
        if let Some(h) = atom.explicit_h {
            return h;
        }

        // Bracket atoms without H spec → 0 (e.g., [Fe], [C])
        // We detect bracket vs organic by whether explicit_h was ever set.
        // But we can't distinguish here without more context stored on the atom.
        // Convention: if element is NOT in organic subset, return 0.
        if !atom.element.is_organic_subset() {
            return 0;
        }

        // Sum bond orders for the atom
        let bond_sum: u8 = self
            .edges(idx)
            .map(|e| e.weight().valence_contribution())
            .sum();

        let valences = atom.element.default_valences();
        if valences.is_empty() || (valences.len() == 1 && valences[0] == 0) {
            return 0;
        }

        if atom.aromatic {
            // For aromatic atoms, the pi system contributes 1 to the effective bond sum.
            // However, if the atom's normal valence is already exactly satisfied by its bonds
            // (e.g., aromatic N bonded to methyl + 2 ring atoms = 3 bonds, valence 3),
            // do NOT add the +1 or we'd incorrectly jump to a higher valence (e.g., N valence 5).
            for &v in valences {
                if v >= bond_sum {
                    let h_raw = v.saturating_sub(bond_sum);
                    let charge_adj = atom.charge.unsigned_abs();
                    if h_raw == 0 {
                        // All valence bonds used; no pi-electron adjustment needed.
                        return 0_u8.saturating_sub(charge_adj);
                    }
                    // Atom still has available valence: apply the +1 aromatic pi adjustment.
                    let effective_sum = bond_sum + 1;
                    for &v2 in valences {
                        if v2 >= effective_sum {
                            return v2.saturating_sub(effective_sum).saturating_sub(charge_adj);
                        }
                    }
                    return 0;
                }
            }
            return 0;
        }

        // Non-aromatic: find smallest valence >= bond_sum
        for &v in valences {
            if v >= bond_sum {
                let charge_adj = atom.charge.unsigned_abs();
                return v.saturating_sub(bond_sum).saturating_sub(charge_adj);
            }
        }

        0
    }

    fn assign_implicit_hydrogens(&mut self) {
        let indices: Vec<NodeIndex> = self.node_indices().collect();
        for idx in indices {
            if self[idx].element == super::atom::Element::H {
                continue; // skip H/D atoms themselves
            }
            let h_implicit = self.compute_implicit_h(idx);
            // Count explicit H/D neighbor atoms in the graph (e.g. [2H] deuterium)
            let h_explicit: u8 = self
                .neighbors(idx)
                .filter(|&nb| self[nb].element == super::atom::Element::H)
                .count()
                .min(255) as u8;
            self[idx].h_count = h_implicit + h_explicit;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::atom::Element;
    use super::*;

    fn make_ethanol() -> MolGraph {
        let mut g = MolGraph::new_undirected();
        let c1 = g.add_node(Atom::new(Element::C));
        let c2 = g.add_node(Atom::new(Element::C));
        let o = g.add_node(Atom::new(Element::O));
        g.add_edge(c1, c2, BondType::Single);
        g.add_edge(c2, o, BondType::Single);
        g.assign_implicit_hydrogens();
        g
    }

    #[test]
    fn test_ethanol_structure() {
        let g = make_ethanol();
        assert_eq!(g.num_atoms(), 3);
        assert_eq!(g.num_bonds(), 2);
    }

    #[test]
    fn test_ethanol_implicit_h() {
        let g = make_ethanol();
        // C-C-O: C1 has 1 heavy neighbor → 3H, C2 has 2 → 2H, O has 1 → 1H
        assert_eq!(g[NodeIndex::new(0)].h_count, 3);
        assert_eq!(g[NodeIndex::new(1)].h_count, 2);
        assert_eq!(g[NodeIndex::new(2)].h_count, 1);
    }

    #[test]
    fn test_benzene_implicit_h() {
        let mut g = MolGraph::new_undirected();
        let atoms: Vec<NodeIndex> = (0..6)
            .map(|_| {
                let mut a = Atom::new(Element::C);
                a.aromatic = true;
                g.add_node(a)
            })
            .collect();
        for i in 0..5 {
            g.add_edge(atoms[i], atoms[i + 1], BondType::Aromatic);
        }
        g.add_edge(atoms[5], atoms[0], BondType::Aromatic);
        g.assign_implicit_hydrogens();
        // Each aromatic C: bond_sum=2 (two aromatic bonds, each contrib 1),
        // effective_sum = 2+1 = 3, valence 4, h = 4-3 = 1
        for &a in &atoms {
            assert_eq!(g[a].h_count, 1, "aromatic C should have 1H");
        }
    }

    #[test]
    fn test_heavy_neighbors() {
        let g = make_ethanol();
        let c2 = NodeIndex::new(1);
        let neighbors = g.heavy_neighbors(c2);
        assert_eq!(neighbors.len(), 2);
    }
}
