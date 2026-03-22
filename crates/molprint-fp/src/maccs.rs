use crate::bitvec::FingerprintBits;
use crate::traits::Fingerprinter;
use molprint_core::{BondType, Element, MolGraph};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

/// Find any simple path of exactly `length` bonds from `start` where the endpoint
/// satisfies `target_pred`. Uses iterative DFS with a per-path visited bitset.
fn has_simple_path(
    mol: &MolGraph,
    start: NodeIndex,
    length: usize,
    target_pred: &dyn Fn(NodeIndex) -> bool,
) -> bool {
    let n = mol.node_count();
    if n <= 64 {
        // Fast path: use single u64 bitmask for visited tracking
        let mut stack: Vec<(NodeIndex, usize, u64)> = Vec::new();
        let init_mask = 1u64 << start.index();
        stack.push((start, length, init_mask));
        while let Some((cur, rem, visited_mask)) = stack.pop() {
            if rem == 0 {
                if target_pred(cur) {
                    return true;
                }
                continue;
            }
            for nb in mol.neighbors(cur) {
                let bit = 1u64 << nb.index();
                if visited_mask & bit == 0 {
                    stack.push((nb, rem - 1, visited_mask | bit));
                }
            }
        }
        false
    } else {
        // Large molecule: use Vec<u64> bitmask per stack frame for correct injectivity
        let nwords = n.div_ceil(64);
        let mut init_mask = vec![0u64; nwords];
        init_mask[start.index() / 64] |= 1u64 << (start.index() % 64);
        let mut stack: Vec<(NodeIndex, usize, Vec<u64>)> = vec![(start, length, init_mask)];
        while let Some((cur, rem, mask)) = stack.pop() {
            if rem == 0 {
                if target_pred(cur) {
                    return true;
                }
                continue;
            }
            for nb in mol.neighbors(cur) {
                let word = nb.index() / 64;
                let bit = 1u64 << (nb.index() % 64);
                if mask[word] & bit == 0 {
                    let mut new_mask = mask.clone();
                    new_mask[word] |= bit;
                    stack.push((nb, rem - 1, new_mask));
                }
            }
        }
        false
    }
}

/// MACCS 166 structural keys fingerprint.
///
/// Bit numbering matches RDKit's `GetMACCSKeysFingerprint` exactly (1-indexed, bits 1-166).
/// Bit 0 and bit 166 (`?`) are always 0. Bits are stored at index N in the bitvector.
pub struct Maccs166;

impl Maccs166 {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Maccs166 {
    fn default() -> Self {
        Self::new()
    }
}

impl Fingerprinter for Maccs166 {
    fn name(&self) -> &str {
        "MACCS166"
    }

    fn nbits(&self) -> usize {
        166
    }

    fn fingerprint(&self, mol: &MolGraph) -> FingerprintBits {
        let mut fp = FingerprintBits::new(166);

        // ── Atom predicates ──────────────────────────────────────────────────────
        // is_heavy: non-hydrogen atom
        let is_heavy = |n: petgraph::graph::NodeIndex| mol[n].element != Element::H;
        // is_heteroatom: non-C, non-H
        let is_hetero = |n: petgraph::graph::NodeIndex| {
            mol[n].element != Element::H && mol[n].element != Element::C
        };
        // "heteroatom with H" = [!#6;!#1;!H0]
        let is_hetero_with_h = |n: petgraph::graph::NodeIndex| {
            mol[n].element != Element::H && mol[n].element != Element::C && mol[n].h_count > 0
        };
        let is_carbon = |n: petgraph::graph::NodeIndex| mol[n].element == Element::C;
        let is_nitrogen = |n: petgraph::graph::NodeIndex| mol[n].element == Element::N;
        let is_oxygen = |n: petgraph::graph::NodeIndex| mol[n].element == Element::O;
        let is_sulfur = |n: petgraph::graph::NodeIndex| mol[n].element == Element::S;
        let is_halogen = |n: petgraph::graph::NodeIndex| {
            matches!(
                mol[n].element,
                Element::F | Element::Cl | Element::Br | Element::I
            )
        };
        let is_aromatic_atom = |n: petgraph::graph::NodeIndex| mol[n].aromatic;
        let in_ring = |n: petgraph::graph::NodeIndex| mol[n].in_ring;
        let h_count = |n: petgraph::graph::NodeIndex| mol[n].h_count;

        // CH2 = carbon with exactly 2 H
        let is_ch2 =
            |n: petgraph::graph::NodeIndex| mol[n].element == Element::C && mol[n].h_count == 2;
        // CH3 = carbon with 3 H
        let is_ch3 =
            |n: petgraph::graph::NodeIndex| mol[n].element == Element::C && mol[n].h_count == 3;
        // [C;H3,H4] = methyl/methane (terminal carbon)
        let is_ch3_or_ch4 = |n: petgraph::graph::NodeIndex| {
            mol[n].element == Element::C && (mol[n].h_count == 3 || mol[n].h_count == 4)
        };
        // NH2 = nitrogen with 2 H
        let is_nh2 =
            |n: petgraph::graph::NodeIndex| mol[n].element == Element::N && mol[n].h_count == 2;
        // N with H
        let is_n_with_h =
            |n: petgraph::graph::NodeIndex| mol[n].element == Element::N && mol[n].h_count > 0;
        // O with H = [O;!H0]
        let is_oh =
            |n: petgraph::graph::NodeIndex| mol[n].element == Element::O && mol[n].h_count > 0;

        // ── Precomputed ──────────────────────────────────────────────────────────
        let nodes: Vec<_> = mol.node_indices().collect();
        let edges: Vec<_> = mol.edge_indices().collect();

        // Bond predicates
        let bond_is_single = |e: petgraph::graph::EdgeIndex| mol[e] == BondType::Single;
        let bond_is_double = |e: petgraph::graph::EdgeIndex| mol[e] == BondType::Double;
        let bond_is_triple = |e: petgraph::graph::EdgeIndex| mol[e] == BondType::Triple;
        // aromatic bond (`:`)
        let bond_is_aromatic = |e: petgraph::graph::EdgeIndex| mol[e] == BondType::Aromatic;

        // bond_in_ring: a bond is a ring bond if removing it keeps its endpoints connected.
        // Precompute for all edges via BFS (O(E*(V+E)) but fine for small molecules).
        let ring_bond_set: std::collections::HashSet<petgraph::graph::EdgeIndex> = edges
            .iter()
            .filter(|&&e| {
                let (u, v) = mol.edge_endpoints(e).unwrap();
                if !mol[u].in_ring || !mol[v].in_ring {
                    return false;
                }
                // BFS from u to v excluding edge e
                let mut visited = std::collections::HashSet::new();
                let mut queue = std::collections::VecDeque::new();
                visited.insert(u);
                queue.push_back(u);
                while let Some(cur) = queue.pop_front() {
                    for er in mol.edges(cur) {
                        if er.id() == e {
                            continue;
                        }
                        let nb = if er.source() == cur {
                            er.target()
                        } else {
                            er.source()
                        };
                        if nb == v {
                            return true;
                        }
                        if visited.insert(nb) {
                            queue.push_back(nb);
                        }
                    }
                }
                false
            })
            .copied()
            .collect();
        let bond_in_ring = |e: petgraph::graph::EdgeIndex| ring_bond_set.contains(&e);
        // `~` = any bond (always true)

        let n_oxygen = nodes.iter().filter(|&&n| is_oxygen(n)).count();
        let n_nitrogen = nodes.iter().filter(|&&n| is_nitrogen(n)).count();
        let n_ch3 = nodes.iter().filter(|&&n| is_ch3(n)).count();
        let n_ch3_or_ch4 = nodes.iter().filter(|&&n| is_ch3_or_ch4(n)).count();
        let n_in_ring = nodes.iter().filter(|&&n| in_ring(n)).count();
        let n_hetero_in_ring = nodes
            .iter()
            .filter(|&&n| is_hetero(n) && in_ring(n))
            .count();
        let n_hetero_with_h = nodes.iter().filter(|&&n| is_hetero_with_h(n)).count();

        // Ring sizes
        let has_ring_size = |sz: u8| nodes.iter().any(|&n| mol[n].ring_sizes.contains(&sz));

        // Count occurrences of a pattern (for min_count logic)
        // count_pattern_nodes(pred) = number of atoms matching pred
        let _count_atoms = |pred: &dyn Fn(petgraph::graph::NodeIndex) -> bool| -> usize {
            nodes.iter().filter(|&&n| pred(n)).count()
        };

        // Any bond between two atom predicates (any bond type)
        let has_bond = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                        pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> bool {
            edges.iter().any(|&e| {
                let (u, v) = mol.edge_endpoints(e).unwrap();
                (pa(u) && pb(v)) || (pa(v) && pb(u))
            })
        };

        // Single bond between atom predicates
        let has_single_bond = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                               pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> bool {
            edges.iter().any(|&e| {
                if !bond_is_single(e) {
                    return false;
                }
                let (u, v) = mol.edge_endpoints(e).unwrap();
                (pa(u) && pb(v)) || (pa(v) && pb(u))
            })
        };

        // Double bond between atom predicates
        let has_double_bond = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                               pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> bool {
            edges.iter().any(|&e| {
                if !bond_is_double(e) {
                    return false;
                }
                let (u, v) = mol.edge_endpoints(e).unwrap();
                (pa(u) && pb(v)) || (pa(v) && pb(u))
            })
        };

        // Triple bond between atom predicates
        let has_triple_bond = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                               pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> bool {
            edges.iter().any(|&e| {
                if !bond_is_triple(e) {
                    return false;
                }
                let (u, v) = mol.edge_endpoints(e).unwrap();
                (pa(u) && pb(v)) || (pa(v) && pb(u))
            })
        };

        // Path A~B~C: find center atom B matching pb, with neighbor matching pa AND neighbor matching pc
        // (pa and pc can be the same atom, we require two distinct neighbors)
        let has_path3 = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                         pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                         pc: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> bool {
            nodes.iter().filter(|&&b| pb(b)).any(|&b| {
                let nbrs: Vec<_> = mol.neighbors(b).collect();
                let a_matches: Vec<_> = nbrs.iter().filter(|&&n| pa(n)).collect();
                let c_matches: Vec<_> = nbrs.iter().filter(|&&n| pc(n)).collect();
                if a_matches.is_empty() || c_matches.is_empty() {
                    return false;
                }
                // Need at least one pair (a,c) with a != c (if pa==pc predicates could overlap)
                a_matches
                    .iter()
                    .any(|&&a| c_matches.iter().any(|&&c| a != c))
            })
        };

        // Path A~B~C where B has the constraint (B center); also handles A~B(~C)~D (3 branches from B)
        let has_path3_branched = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                                  pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                                  pc: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                                  pd: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> bool {
            // A~B(~C)~D: B has neighbors matching A, C, D (all distinct)
            nodes.iter().filter(|&&b| pb(b)).any(|&b| {
                let nbrs: Vec<_> = mol.neighbors(b).collect();
                nbrs.iter().any(|&a| {
                    if !pa(a) {
                        return false;
                    }
                    nbrs.iter().any(|&c| {
                        if c == a || !pc(c) {
                            return false;
                        }
                        nbrs.iter().any(|&d| d != a && d != c && pd(d))
                    })
                })
            })
        };

        // Path A~B~C~D (4 atoms, linear)
        let has_path4 = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                         pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                         pc: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                         pd: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> bool {
            edges.iter().any(|&e| {
                let (b, c) = mol.edge_endpoints(e).unwrap();
                // Try b=B, c=C; require A≠D for full injectivity
                if pb(b) && pc(c) {
                    for a in mol.neighbors(b) {
                        if a == c || !pa(a) {
                            continue;
                        }
                        for d in mol.neighbors(c) {
                            if d == b || d == a || !pd(d) {
                                continue;
                            }
                            return true;
                        }
                    }
                }
                // Try b=C, c=B (reverse)
                if pb(c) && pc(b) {
                    for a in mol.neighbors(c) {
                        if a == b || !pa(a) {
                            continue;
                        }
                        for d in mol.neighbors(b) {
                            if d == c || d == a || !pd(d) {
                                continue;
                            }
                            return true;
                        }
                    }
                }
                false
            })
        };

        // Path A~*~*~*~B (5 atoms)
        // has_path5: find any SIMPLE path of exactly 4 bonds from an A-atom to a B-atom.
        // Uses DFS (not BFS) so it correctly finds non-shortest paths in rings.
        let has_path5 = |pa: &dyn Fn(NodeIndex) -> bool, pe: &dyn Fn(NodeIndex) -> bool| -> bool {
            nodes
                .iter()
                .filter(|&&n| pa(n))
                .any(|&start| has_simple_path(mol, start, 4, pe))
        };

        // Count occurrences of A~B~* patterns (edge patterns for min_count)
        // counts how many atoms satisfy: is_pred(n) and has neighbor satisfying other_pred
        let count_bonded = |pa: &dyn Fn(petgraph::graph::NodeIndex) -> bool,
                            pb: &dyn Fn(petgraph::graph::NodeIndex) -> bool|
         -> usize {
            // Count unique edges satisfying pa~pb
            edges
                .iter()
                .filter(|&&e| {
                    let (u, v) = mol.edge_endpoints(e).unwrap();
                    (pa(u) && pb(v)) || (pa(v) && pb(u))
                })
                .count()
        };

        // ── Keys 1-166 ──────────────────────────────────────────────────────────
        // Key 1: `?` always 0
        // Key 2: [#104] (Rutherfordium) — not in our element set, always 0

        // Key 3: [#32,#33,#34,#50,#51,#52,#82,#83,#84] = Ge,As,Se,Sn,Sb,Te,Pb,Bi,Po
        if nodes.iter().any(|&n| {
            matches!(
                mol[n].element,
                Element::Ge
                    | Element::As
                    | Element::Se
                    | Element::Sn
                    | Element::Sb
                    | Element::Te
                    | Element::Pb
                    | Element::Bi
            )
        }) {
            fp.set(3);
        }

        // Key 4: [Ac,Th,Pa,U,...] actinides — not in our element set, always 0

        // Key 5: [Sc,Ti,Y,Zr,Hf]
        if nodes
            .iter()
            .any(|&n| matches!(mol[n].element, Element::Ti | Element::Zr))
        {
            fp.set(5);
        }

        // Key 6: lanthanides — not in our element set, always 0

        // Key 7: [V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]
        if nodes.iter().any(|&n| {
            matches!(
                mol[n].element,
                Element::V | Element::Cr | Element::Mn | Element::Mo
            )
        }) {
            fp.set(7);
        }

        // Key 8: [!#6;!#1]1~*~*~*~1 = 4-membered ring containing a heteroatom
        if has_ring_size(4)
            && nodes
                .iter()
                .any(|&n| is_hetero(n) && in_ring(n) && mol[n].ring_sizes.contains(&4))
        {
            fp.set(8);
        }

        // Key 9: [Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]
        if nodes.iter().any(|&n| {
            matches!(
                mol[n].element,
                Element::Fe
                    | Element::Co
                    | Element::Ni
                    | Element::Ru
                    | Element::Rh
                    | Element::Pd
                    | Element::Pt
            )
        }) {
            fp.set(9);
        }

        // Key 10: [Be,Mg,Ca,Sr,Ba,Ra]
        if nodes.iter().any(|&n| {
            matches!(
                mol[n].element,
                Element::Mg | Element::Ca | Element::Sr | Element::Ba
            )
        }) {
            fp.set(10);
        }

        // Key 11: *1~*~*~*~1 = any 4-membered ring
        if has_ring_size(4) {
            fp.set(11);
        }

        // Key 12: [Cu,Zn,Ag,Cd,Au,Hg]
        if nodes.iter().any(|&n| {
            matches!(
                mol[n].element,
                Element::Cu | Element::Zn | Element::Ag | Element::Cd | Element::Au | Element::Hg
            )
        }) {
            fp.set(12);
        }

        // Key 13: [#8]~[#7](~[#6])~[#6] = O-N(-C)-C
        if has_path3_branched(&is_oxygen, &is_nitrogen, &is_carbon, &is_carbon) {
            fp.set(13);
        }

        // Key 14: [#16]-[#16] = S-S single bond
        if has_single_bond(&is_sulfur, &is_sulfur) {
            fp.set(14);
        }

        // Key 15: [#8]~[#6](~[#8])~[#8] = C bonded to 3 O (carbonate)
        if nodes
            .iter()
            .filter(|&&n| is_carbon(n))
            .any(|&c| mol.neighbors(c).filter(|&nb| is_oxygen(nb)).count() >= 3)
        {
            fp.set(15);
        }

        // Key 16: [!#6;!#1]1~*~*~1 = 3-membered ring with heteroatom
        if has_ring_size(3)
            && nodes
                .iter()
                .any(|&n| is_hetero(n) && in_ring(n) && mol[n].ring_sizes.contains(&3))
        {
            fp.set(16);
        }

        // Key 17: [#6]#[#6] = C≡C
        if has_triple_bond(&is_carbon, &is_carbon) {
            fp.set(17);
        }

        // Key 18: [#5,#13,#31,#49,#81] = B,Al,Ga,In,Tl
        if nodes.iter().any(|&n| {
            matches!(
                mol[n].element,
                Element::B | Element::Al | Element::Ga | Element::In | Element::Tl
            )
        }) {
            fp.set(18);
        }

        // Key 19: *1~*~*~*~*~*~*~1 = 7-membered ring (includes non-SSSR outer rings)
        {
            let has_7_ring = has_ring_size(7) || {
                let n = mol.node_count();
                let mut found = false;
                'dfs19: for &start in &nodes {
                    if !in_ring(start) {
                        continue;
                    }
                    if n <= 64 {
                        let init_vis: u64 = 1u64 << start.index();
                        let mut stack: Vec<(NodeIndex, usize, u64)> = Vec::new();
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            stack.push((nb, 1, init_vis | (1u64 << nb.index())));
                        }
                        while let Some((cur, depth, vis)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 6 {
                                        found = true;
                                        break 'dfs19;
                                    }
                                    continue;
                                }
                                let bit = 1u64 << nb.index();
                                if vis & bit != 0 {
                                    continue;
                                }
                                if depth < 6 {
                                    stack.push((nb, depth + 1, vis | bit));
                                }
                            }
                        }
                    } else {
                        let mut stack: Vec<(
                            NodeIndex,
                            usize,
                            std::collections::HashSet<NodeIndex>,
                        )> = Vec::new();
                        let mut init_vis = std::collections::HashSet::new();
                        init_vis.insert(start);
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            let mut vis = init_vis.clone();
                            vis.insert(nb);
                            stack.push((nb, 1, vis));
                        }
                        while let Some((cur, depth, vis)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 6 {
                                        found = true;
                                        break 'dfs19;
                                    }
                                    continue;
                                }
                                if vis.contains(&nb) {
                                    continue;
                                }
                                if depth < 6 {
                                    let mut new_vis = vis.clone();
                                    new_vis.insert(nb);
                                    stack.push((nb, depth + 1, new_vis));
                                }
                            }
                        }
                    }
                }
                found
            };
            if has_7_ring {
                fp.set(19);
            }
        }

        // Key 20: [#14] = Si
        if nodes.iter().any(|&n| mol[n].element == Element::Si) {
            fp.set(20);
        }

        // Key 21: [#6]=[#6](~[!#6;!#1])~[!#6;!#1] = C=C where ONE C has ≥2 heteroatom neighbors
        // Both endpoints must be carbon; one of them has ≥2 heteroatom neighbors besides the other.
        if edges.iter().any(|&e| {
            if !bond_is_double(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            if !is_carbon(u) || !is_carbon(v) {
                return false;
            }
            let u_het = mol
                .neighbors(u)
                .filter(|&nb| nb != v && is_hetero(nb))
                .count();
            let v_het = mol
                .neighbors(v)
                .filter(|&nb| nb != u && is_hetero(nb))
                .count();
            u_het >= 2 || v_het >= 2
        }) {
            fp.set(21);
        }

        // Key 22: *1~*~*~1 = 3-membered ring
        if has_ring_size(3) {
            fp.set(22);
        }

        // Key 23: [#7]~[#6](~[#8])~[#8] = N-C(-O)-O
        if has_path3_branched(&is_nitrogen, &is_carbon, &is_oxygen, &is_oxygen) {
            fp.set(23);
        }

        // Key 24: [#7]-[#8] = N-O single bond
        if has_single_bond(&is_nitrogen, &is_oxygen) {
            fp.set(24);
        }

        // Key 25: [#7]~[#6](~[#7])~[#7] = N-C(-N)-N (guanidine-like)
        if has_path3_branched(&is_nitrogen, &is_carbon, &is_nitrogen, &is_nitrogen) {
            fp.set(25);
        }

        // Key 26: [#6]=;@[#6](@*)@* = C=C ring bond, second C has 2 more ring bonds
        // The C=C bond itself must be a ring bond (=;@), and the second C must have
        // 2 additional ring bonds (@*)@*  — i.e. the second C is at a ring junction.
        if edges.iter().any(|&e| {
            if !bond_is_double(e) || !bond_in_ring(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            if !is_carbon(u) || !is_carbon(v) {
                return false;
            }
            // Check if either endpoint has 2+ additional ring bonds (ring junction)
            let extra_ring_bonds = |center: NodeIndex, other: NodeIndex| -> usize {
                mol.edges(center)
                    .filter(|er| er.id() != e && ring_bond_set.contains(&er.id()))
                    .filter(|er| {
                        let nb = if er.source() == center {
                            er.target()
                        } else {
                            er.source()
                        };
                        nb != other
                    })
                    .count()
            };
            extra_ring_bonds(u, v) >= 2 || extra_ring_bonds(v, u) >= 2
        }) {
            fp.set(26);
        }

        // Key 27: [I] = Iodine
        if nodes.iter().any(|&n| mol[n].element == Element::I) {
            fp.set(27);
        }

        // Key 28: [!#6;!#1]~[CH2]~[!#6;!#1] = heteroatom-CH2-heteroatom
        if has_path3(&is_hetero, &is_ch2, &is_hetero) {
            fp.set(28);
        }

        // Key 29: [#15] = P
        if nodes.iter().any(|&n| mol[n].element == Element::P) {
            fp.set(29);
        }

        // Key 30: [#6]~[!#6;!#1](~[#6])(~[#6])~* = heteroatom with 3 C neighbors AND 1 more (degree>=4)
        if nodes.iter().filter(|&&n| is_hetero(n)).any(|&n| {
            mol.neighbors(n).filter(|&nb| is_carbon(nb)).count() >= 3
                && mol.neighbors(n).count() >= 4
        }) {
            fp.set(30);
        }

        // Key 31: [!#6;!#1]~[F,Cl,Br,I] = heteroatom bonded to halogen
        if has_bond(&is_hetero, &is_halogen) {
            fp.set(31);
        }

        // Key 32: [#6]~[#16]~[#7] = C-S-N
        if has_path3(&is_carbon, &is_sulfur, &is_nitrogen) {
            fp.set(32);
        }

        // Key 33: [#7]~[#16] = N-S bond (any)
        if has_bond(&is_nitrogen, &is_sulfur) {
            fp.set(33);
        }

        // Key 34: [CH2]=* = CH2 in a double bond
        if edges.iter().any(|&e| {
            bond_is_double(e) && {
                let (u, v) = mol.edge_endpoints(e).unwrap();
                is_ch2(u) || is_ch2(v)
            }
        }) {
            fp.set(34);
        }

        // Key 35: [Li,Na,K,Rb,Cs,Fr] = alkali metals
        if nodes.iter().any(|&n| {
            matches!(
                mol[n].element,
                Element::Li | Element::Na | Element::K | Element::Rb | Element::Cs
            )
        }) {
            fp.set(35);
        }

        // Key 36: [#16R] = S in ring
        if nodes.iter().any(|&n| is_sulfur(n) && in_ring(n)) {
            fp.set(36);
        }

        // Key 37: [#7]~[#6](~[#8])~[#7] = N-C(-O)-N
        if has_path3_branched(&is_nitrogen, &is_carbon, &is_oxygen, &is_nitrogen) {
            fp.set(37);
        }

        // Key 38: [#7]~[#6](~[#6])~[#7] = N-C(-C)-N
        if has_path3_branched(&is_nitrogen, &is_carbon, &is_carbon, &is_nitrogen) {
            fp.set(38);
        }

        // Key 39: [#8]~[#16](~[#8])~[#8] = O-S(-O)-O (sulfate-like)
        if nodes
            .iter()
            .filter(|&&n| is_sulfur(n))
            .any(|&s| mol.neighbors(s).filter(|&nb| is_oxygen(nb)).count() >= 3)
        {
            fp.set(39);
        }

        // Key 40: [#16]-[#8] = S-O single bond
        if has_single_bond(&is_sulfur, &is_oxygen) {
            fp.set(40);
        }

        // Key 41: [#6]#[#7] = C≡N (nitrile)
        if has_triple_bond(&is_carbon, &is_nitrogen) {
            fp.set(41);
        }

        // Key 42: F = Fluorine
        if nodes.iter().any(|&n| mol[n].element == Element::F) {
            fp.set(42);
        }

        // Key 43: [!#6;!#1;!H0]~*~[!#6;!#1;!H0] = two H-bearing heteroatoms 1 apart
        if has_path3(&is_hetero_with_h, &is_heavy, &is_hetero_with_h) {
            fp.set(43);
        }

        // Key 44: [!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53] = uncommon heteroatom
        // Not H(1), C(6), N(7), O(8), F(9), Si(14), P(15), S(16), Cl(17), Br(35), I(53)
        if nodes.iter().any(|&n| {
            !matches!(
                mol[n].element,
                Element::H
                    | Element::C
                    | Element::N
                    | Element::O
                    | Element::F
                    | Element::Si
                    | Element::P
                    | Element::S
                    | Element::Cl
                    | Element::Br
                    | Element::I
            )
        }) {
            fp.set(44);
        }

        // Key 45: [#6]=[#6]~[#7] = C=C-N
        if has_path3(&is_carbon, &is_carbon, &is_nitrogen)
            && has_double_bond(&is_carbon, &is_carbon)
        {
            // More precise: find C=C where one C also has N neighbor
            let found = edges.iter().any(|&e| {
                if !bond_is_double(e) {
                    return false;
                }
                let (u, v) = mol.edge_endpoints(e).unwrap();
                if !is_carbon(u) || !is_carbon(v) {
                    return false;
                }
                mol.neighbors(u).any(|nb| nb != v && is_nitrogen(nb))
                    || mol.neighbors(v).any(|nb| nb != u && is_nitrogen(nb))
            });
            if found {
                fp.set(45);
            }
        }

        // Key 46: Br = Bromine
        if nodes.iter().any(|&n| mol[n].element == Element::Br) {
            fp.set(46);
        }

        // Key 47: [#16]~*~[#7] = S-*-N
        if has_path3(&is_sulfur, &is_heavy, &is_nitrogen) {
            fp.set(47);
        }

        // Key 48: [#8]~[!#6;!#1](~[#8])(~[#8]) = O-X(-O)-O (3 oxygens on heteroatom)
        if nodes
            .iter()
            .filter(|&&n| is_hetero(n))
            .any(|&x| mol.neighbors(x).filter(|&nb| is_oxygen(nb)).count() >= 3)
        {
            fp.set(48);
        }

        // Key 49: [!+0] = charged atom
        if nodes.iter().any(|&n| mol[n].charge != 0) {
            fp.set(49);
        }

        // Key 50: [#6]=[#6](~[#6])~[#6] = C=C where one C has ≥2 extra C neighbors
        // SMARTS: second C has a C branch AND a C continuation = 2 extra C neighbors
        if edges.iter().any(|&e| {
            if !bond_is_double(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            if !is_carbon(u) || !is_carbon(v) {
                return false;
            }
            let u_extra_c = mol
                .neighbors(u)
                .filter(|&nb| nb != v && is_carbon(nb))
                .count();
            let v_extra_c = mol
                .neighbors(v)
                .filter(|&nb| nb != u && is_carbon(nb))
                .count();
            u_extra_c >= 2 || v_extra_c >= 2
        }) {
            fp.set(50);
        }

        // Key 51: [#6]~[#16]~[#8] = C-S-O
        if has_path3(&is_carbon, &is_sulfur, &is_oxygen) {
            fp.set(51);
        }

        // Key 52: [#7]~[#7] = N-N
        if has_bond(&is_nitrogen, &is_nitrogen) {
            fp.set(52);
        }

        // Key 53: [!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0] = two H-bearing heteroatoms 3 bonds apart
        if has_path5(&is_hetero_with_h, &is_hetero_with_h) {
            fp.set(53);
        }

        // Key 54: [!#6;!#1;!H0]~*~*~[!#6;!#1;!H0] = two H-bearing heteroatoms 2 bonds apart
        if has_path4(&is_hetero_with_h, &is_heavy, &is_heavy, &is_hetero_with_h) {
            fp.set(54);
        }

        // Key 55: [#8]~[#16]~[#8] = O-S-O
        if has_path3(&is_oxygen, &is_sulfur, &is_oxygen) {
            fp.set(55);
        }

        // Key 56: [#8]~[#7](~[#8])~[#6] = O-N(-O)-C (nitro group)
        if has_path3_branched(&is_oxygen, &is_nitrogen, &is_oxygen, &is_carbon) {
            fp.set(56);
        }

        // Key 57: [#8R] = O in ring
        if nodes.iter().any(|&n| is_oxygen(n) && in_ring(n)) {
            fp.set(57);
        }

        // Key 58: [!#6;!#1]~[#16]~[!#6;!#1] = heteroatom-S-heteroatom
        if has_path3(&is_hetero, &is_sulfur, &is_hetero) {
            fp.set(58);
        }

        // Key 59: [#16]!:*:* = S-X where X is in an aromatic ring (non-aromatic bond S to aromatic atom)
        if edges.iter().any(|&e| {
            if bond_is_aromatic(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            let other = if is_sulfur(u) {
                v
            } else if is_sulfur(v) {
                u
            } else {
                return false;
            };
            // 'other' must have an aromatic bond (i.e., be in an aromatic ring)
            mol.edges(other).any(|er| bond_is_aromatic(er.id()))
        }) {
            fp.set(59);
        }

        // Key 60: [#16]=[#8] = S=O
        if has_double_bond(&is_sulfur, &is_oxygen) {
            fp.set(60);
        }

        // Key 61: *~[#16](~*)~* = S with 3+ bonds (trisubstituted S)
        if nodes
            .iter()
            .filter(|&&n| is_sulfur(n))
            .any(|&s| mol.neighbors(s).count() >= 3)
        {
            fp.set(61);
        }

        // Key 62: *@*!@*@* = two ring atoms bridged by a non-ring bond (ring fusion / spiro pattern)
        // Approximate: any atom in a ring that has a non-ring neighbor that is also in a ring
        if edges.iter().any(|&e| {
            if bond_in_ring(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            in_ring(u) && in_ring(v)
        }) {
            fp.set(62);
        }

        // Key 63: [#7]=[#8] = N=O
        if has_double_bond(&is_nitrogen, &is_oxygen) {
            fp.set(63);
        }

        // Key 64: *@*!@[#16] = ring atom bonded to S via non-ring bond
        if edges.iter().any(|&e| {
            if bond_in_ring(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            (in_ring(u) && is_sulfur(v)) || (in_ring(v) && is_sulfur(u))
        }) {
            fp.set(64);
        }

        // Key 65: c:n = aromatic C bonded to aromatic N via aromatic bond
        if edges.iter().any(|&e| {
            bond_is_aromatic(e) && {
                let (u, v) = mol.edge_endpoints(e).unwrap();
                (is_carbon(u) && is_nitrogen(v)) || (is_carbon(v) && is_nitrogen(u))
            }
        }) {
            fp.set(65);
        }

        // Key 66: [#6]~[#6](~[#6])(~[#6])~* = C bonded to 3+ C neighbors AND degree>=4
        if nodes.iter().filter(|&&n| is_carbon(n)).any(|&c| {
            mol.neighbors(c).filter(|&nb| is_carbon(nb)).count() >= 3
                && mol.neighbors(c).count() >= 4
        }) {
            fp.set(66);
        }

        // Key 67: [!#6;!#1]~[#16] = heteroatom bonded to S
        if has_bond(&is_hetero, &is_sulfur) {
            fp.set(67);
        }

        // Key 68: [!#6;!#1;!H0]~[!#6;!#1;!H0] = two H-bearing heteroatoms directly bonded
        if has_bond(&is_hetero_with_h, &is_hetero_with_h) {
            fp.set(68);
        }

        // Key 69: [!#6;!#1]~[!#6;!#1;!H0] = heteroatom bonded to H-bearing heteroatom
        if has_bond(&is_hetero, &is_hetero_with_h) {
            fp.set(69);
        }

        // Key 70: [!#6;!#1]~[#7]~[!#6;!#1] = heteroatom-N-heteroatom
        if has_path3(&is_hetero, &is_nitrogen, &is_hetero) {
            fp.set(70);
        }

        // Key 71: [#7]~[#8] = N-O (any bond)
        if has_bond(&is_nitrogen, &is_oxygen) {
            fp.set(71);
        }

        // Key 72: [#8]~*~*~[#8] = O-*-*-O (two oxygens 2 bonds apart)
        if has_path4(&is_oxygen, &is_heavy, &is_heavy, &is_oxygen) {
            fp.set(72);
        }

        // Key 73: [#16]=* = S in a double bond
        if edges.iter().any(|&e| {
            bond_is_double(e) && {
                let (u, v) = mol.edge_endpoints(e).unwrap();
                is_sulfur(u) || is_sulfur(v)
            }
        }) {
            fp.set(73);
        }

        // Key 74: [CH3]~*~[CH3] = two CH3 groups one atom apart
        if has_path3(&is_ch3, &is_heavy, &is_ch3) {
            fp.set(74);
        }

        // Key 75: *!@[#7]@* = N in ring bonded to non-ring atom
        // N in ring with at least one non-ring-bond neighbor
        if nodes
            .iter()
            .filter(|&&n| is_nitrogen(n) && in_ring(n))
            .any(|&n| {
                edges.iter().any(|&e| {
                    if bond_in_ring(e) {
                        return false;
                    }
                    let (u, v) = mol.edge_endpoints(e).unwrap();
                    u == n || v == n
                })
            })
        {
            fp.set(75);
        }

        // Key 76: [#6]=[#6](~*)~* = C=C where one C has ≥2 extra heavy-atom neighbors
        // SMARTS: second C has a branch (~*) AND a continuation (~*) = 2 extra neighbors
        if edges.iter().any(|&e| {
            if !bond_is_double(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            if !is_carbon(u) || !is_carbon(v) {
                return false;
            }
            let u_extra = mol.neighbors(u).filter(|&nb| nb != v).count();
            let v_extra = mol.neighbors(v).filter(|&nb| nb != u).count();
            u_extra >= 2 || v_extra >= 2
        }) {
            fp.set(76);
        }

        // Key 77: [#7]~*~[#7] = N-*-N (two N two bonds apart)
        if has_path3(&is_nitrogen, &is_heavy, &is_nitrogen) {
            fp.set(77);
        }

        // Key 78: [#6]=[#7] = C=N
        if has_double_bond(&is_carbon, &is_nitrogen) {
            fp.set(78);
        }

        // Key 79: [#7]~*~*~[#7] = N-*-*-N (two N three bonds apart)
        if has_path4(&is_nitrogen, &is_heavy, &is_heavy, &is_nitrogen) {
            fp.set(79);
        }

        // Key 80: [#7]~*~*~*~[#7] = N-*-*-*-N (two N four bonds apart)
        if has_path5(&is_nitrogen, &is_nitrogen) {
            fp.set(80);
        }

        // Key 81: [#16]~*(~*)~* = S bonded to an atom that has 2+ other heavy-atom neighbors
        // SMARTS: S~center(~*)~*  → center has S + 2 more bonds = degree ≥ 3
        if nodes.iter().filter(|&&n| is_sulfur(n)).any(|&s| {
            mol.neighbors(s).any(|center| {
                // center must have ≥2 neighbors besides s (i.e. degree ≥ 3 total)
                mol.neighbors(center)
                    .filter(|&nb| nb != s && is_heavy(nb))
                    .count()
                    >= 2
            })
        }) {
            fp.set(81);
        }

        // Key 82: *~[CH2]~[!#6;!#1;!H0] = CH2 bonded to H-bearing heteroatom
        if has_path3(&is_heavy, &is_ch2, &is_hetero_with_h) {
            fp.set(82);
        }

        // Key 83: [!#6;!#1]1~*~*~*~*~1 = 5-membered ring with heteroatom
        if has_ring_size(5)
            && nodes
                .iter()
                .any(|&n| is_hetero(n) && in_ring(n) && mol[n].ring_sizes.contains(&5))
        {
            fp.set(83);
        }

        // Key 84: [NH2] = NH2
        if nodes.iter().any(|&n| is_nh2(n)) {
            fp.set(84);
        }

        // Key 85: [#6]~[#7](~[#6])~[#6] = tertiary amine (N with 3 C neighbors)
        if nodes
            .iter()
            .filter(|&&n| is_nitrogen(n))
            .any(|&n| mol.neighbors(n).filter(|&nb| is_carbon(nb)).count() >= 3)
        {
            fp.set(85);
        }

        // Key 86: [C;H2,H3][!#6;!#1][C;H2,H3] = CH2/CH3-heteroatom-CH2/CH3
        {
            let is_ch2_or_ch3 = |n: petgraph::graph::NodeIndex| {
                is_carbon(n) && (h_count(n) == 2 || h_count(n) == 3)
            };
            if has_path3(&is_ch2_or_ch3, &is_hetero, &is_ch2_or_ch3) {
                fp.set(86);
            }
        }

        // Key 87: [F,Cl,Br,I]!@*@* = halogen outside ring bonded to ring atom
        if edges.iter().any(|&e| {
            if bond_in_ring(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            (is_halogen(u) && in_ring(v)) || (is_halogen(v) && in_ring(u))
        }) {
            fp.set(87);
        }

        // Key 88: [#16] = any S
        if nodes.iter().any(|&n| is_sulfur(n)) {
            fp.set(88);
        }

        // Key 89: [#8]~*~*~*~[#8] = O-*-*-*-O (two O four bonds apart)
        if has_path5(&is_oxygen, &is_oxygen) {
            fp.set(89);
        }

        // Key 90: [!#6;!#1;!H0]~*~*~[CH2]~* = 5-atom path; CH2 must have an outer heavy neighbor.
        // has_path4 only finds 4 atoms; the trailing * requires CH2 to have a 2nd heavy bond.
        if edges.iter().any(|&e| {
            let (b, c) = mol.edge_endpoints(e).unwrap();
            let try_dir = |b_node: NodeIndex, c_node: NodeIndex| -> bool {
                if !is_ch2(c_node) {
                    return false;
                }
                // CH2 must have outer neighbor ≠ b_node (the trailing *)
                if !mol.neighbors(c_node).any(|x| x != b_node) {
                    return false;
                }
                // Find A~mid~b_node where A is hetero_with_h and mid ≠ c_node
                mol.neighbors(b_node)
                    .filter(|&mid| mid != c_node)
                    .any(|mid| {
                        mol.neighbors(mid)
                            .any(|a| a != b_node && is_hetero_with_h(a))
                    })
            };
            try_dir(b, c) || try_dir(c, b)
        }) {
            fp.set(90);
        }

        // Key 91: [$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@...)] min_count=0
        // Pattern 1: [hetero_with_h]~*~*~*~[CH2]~* (6 atoms, 5 bonds)
        // The CH2 must have a further neighbor (trailing ~*), so graph degree >= 2.
        // Terminal vinyl CH2 (degree=1 in graph) must be excluded.
        {
            let found91_p1 = nodes.iter().any(|&start| {
                if !is_hetero_with_h(start) {
                    return false;
                }
                has_simple_path(mol, start, 4, &|ch2_node| {
                    if !is_ch2(ch2_node) {
                        return false;
                    }
                    // Trailing ~*: ch2 must have ≥2 graph neighbors (one from path, one more)
                    mol.neighbors(ch2_node).count() >= 2
                })
            });
            if found91_p1 {
                fp.set(91);
            }
        }

        // Key 92: [#8]~[#6](~[#7])~[#6] = O-C(-N)-C (amide-like)
        if has_path3_branched(&is_oxygen, &is_carbon, &is_nitrogen, &is_carbon) {
            fp.set(92);
        }

        // Key 93: [!#6;!#1]~[CH3] = heteroatom bonded to CH3
        if has_bond(&is_hetero, &is_ch3) {
            fp.set(93);
        }

        // Key 94: [!#6;!#1]~[#7] = heteroatom bonded to N
        if has_bond(&is_hetero, &is_nitrogen) {
            fp.set(94);
        }

        // Key 95: [#7]~*~*~[#8] = N-*-*-O (N and O 3 bonds apart)
        if has_path4(&is_nitrogen, &is_heavy, &is_heavy, &is_oxygen) {
            fp.set(95);
        }

        // Key 96: *1~*~*~*~*~1 = any 5-membered ring
        if has_ring_size(5) {
            fp.set(96);
        }

        // Key 97: [#7]~*~*~*~[#8] = N-*-*-*-O (N and O 4 bonds apart)
        if has_path5(&is_nitrogen, &is_oxygen) {
            fp.set(97);
        }

        // Key 98: [!#6;!#1]1~*~*~*~*~*~1 = 6-membered ring with heteroatom (incl. non-SSSR)
        {
            let has98 = (has_ring_size(6)
                && nodes
                    .iter()
                    .any(|&n| is_hetero(n) && in_ring(n) && mol[n].ring_sizes.contains(&6)))
                || {
                    // DFS over ring-bond subgraph: find 6-cycle containing a heteroatom
                    let n = mol.node_count();
                    let mut found = false;
                    'dfs98: for &start in &nodes {
                        if !in_ring(start) || !is_hetero(start) {
                            continue;
                        }
                        if n <= 64 {
                            let init_vis: u64 = 1u64 << start.index();
                            let mut stack: Vec<(NodeIndex, usize, u64)> = Vec::new();
                            for er in mol.edges(start) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == start {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                stack.push((nb, 1, init_vis | (1u64 << nb.index())));
                            }
                            while let Some((cur, depth, vis)) = stack.pop() {
                                for er in mol.edges(cur) {
                                    if !ring_bond_set.contains(&er.id()) {
                                        continue;
                                    }
                                    let nb = if er.source() == cur {
                                        er.target()
                                    } else {
                                        er.source()
                                    };
                                    if nb == start {
                                        if depth == 5 {
                                            found = true;
                                            break 'dfs98;
                                        }
                                        continue;
                                    }
                                    let bit = 1u64 << nb.index();
                                    if vis & bit != 0 {
                                        continue;
                                    }
                                    if depth < 5 {
                                        stack.push((nb, depth + 1, vis | bit));
                                    }
                                }
                            }
                        } else {
                            let mut stack: Vec<(
                                NodeIndex,
                                usize,
                                std::collections::HashSet<NodeIndex>,
                            )> = Vec::new();
                            let mut init_vis = std::collections::HashSet::new();
                            init_vis.insert(start);
                            for er in mol.edges(start) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == start {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                let mut vis = init_vis.clone();
                                vis.insert(nb);
                                stack.push((nb, 1, vis));
                            }
                            while let Some((cur, depth, vis)) = stack.pop() {
                                for er in mol.edges(cur) {
                                    if !ring_bond_set.contains(&er.id()) {
                                        continue;
                                    }
                                    let nb = if er.source() == cur {
                                        er.target()
                                    } else {
                                        er.source()
                                    };
                                    if nb == start {
                                        if depth == 5 {
                                            found = true;
                                            break 'dfs98;
                                        }
                                        continue;
                                    }
                                    if vis.contains(&nb) {
                                        continue;
                                    }
                                    if depth < 5 {
                                        let mut nv = vis.clone();
                                        nv.insert(nb);
                                        stack.push((nb, depth + 1, nv));
                                    }
                                }
                            }
                        }
                    }
                    found
                };
            if has98 {
                fp.set(98);
            }
        }

        // Key 99: [#6]=[#6] = C=C (any)
        if has_double_bond(&is_carbon, &is_carbon) {
            fp.set(99);
        }

        // Key 100: *~[CH2]~[#7] = CH2 bonded to N
        if has_path3(&is_heavy, &is_ch2, &is_nitrogen) {
            fp.set(100);
        }

        // Key 101: atom in an 8–14 membered ring (includes non-SSSR outer rings, e.g. naphthalene's 10-ring).
        // Check SSSR first (fast path), then DFS over ring-bond subgraph for non-SSSR cycles.
        {
            let has_large_ring = (8u8..=14).any(has_ring_size) || {
                // DFS: find cycles of size 8–14 in the ring-bond subgraph.
                // When we're at depth d and see `start` as a neighbor, cycle size = d + 1.
                // So we need d in [7, 13] for cycle sizes [8, 14].
                let n = mol.node_count();
                let mut found = false;
                'dfs_outer: for &start in &nodes {
                    if !in_ring(start) {
                        continue;
                    }
                    if n <= 64 {
                        // Fast path: bitmask visited set
                        let init_vis: u64 = 1u64 << start.index();
                        let mut stack: Vec<(NodeIndex, usize, u64)> = Vec::new();
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            stack.push((nb, 1, init_vis | (1u64 << nb.index())));
                        }
                        while let Some((cur, depth, vis)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if (7..=13).contains(&depth) {
                                        found = true;
                                        break 'dfs_outer;
                                    }
                                    continue;
                                }
                                let bit = 1u64 << nb.index();
                                if vis & bit != 0 {
                                    continue;
                                }
                                if depth < 13 {
                                    stack.push((nb, depth + 1, vis | bit));
                                }
                            }
                        }
                    } else {
                        // Fallback for large molecules: HashSet visited set
                        let mut stack: Vec<(
                            NodeIndex,
                            usize,
                            std::collections::HashSet<NodeIndex>,
                        )> = Vec::new();
                        let mut init_vis = std::collections::HashSet::new();
                        init_vis.insert(start);
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            let mut vis = init_vis.clone();
                            vis.insert(nb);
                            stack.push((nb, 1, vis));
                        }
                        while let Some((cur, depth, vis)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if (7..=13).contains(&depth) {
                                        found = true;
                                        break 'dfs_outer;
                                    }
                                    continue;
                                }
                                if vis.contains(&nb) {
                                    continue;
                                }
                                if depth < 13 {
                                    let mut new_vis = vis.clone();
                                    new_vis.insert(nb);
                                    stack.push((nb, depth + 1, new_vis));
                                }
                            }
                        }
                    }
                }
                found
            };
            if has_large_ring {
                fp.set(101);
            }
        }

        // Key 102: [!#6;!#1]~[#8] = heteroatom bonded to O
        if has_bond(&is_hetero, &is_oxygen) {
            fp.set(102);
        }

        // Key 103: Cl = Chlorine
        if nodes.iter().any(|&n| mol[n].element == Element::Cl) {
            fp.set(103);
        }

        // Key 104: [!#6;!#1;!H0]~*~[CH2]~* = H-bearing heteroatom-*-CH2-*
        if has_path4(&is_hetero_with_h, &is_heavy, &is_ch2, &is_heavy) {
            fp.set(104);
        }

        // Key 105: *@*(@*)@* = ring atom with 3+ ring bonds (ring branch point)
        if nodes.iter().filter(|&&n| in_ring(n)).any(|&n| {
            edges
                .iter()
                .filter(|&&e| {
                    if !bond_in_ring(e) {
                        return false;
                    }
                    let (u, v) = mol.edge_endpoints(e).unwrap();
                    u == n || v == n
                })
                .count()
                >= 3
        }) {
            fp.set(105);
        }

        // Key 106: [!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1] = atom bonded to 3 heteroatoms
        if nodes
            .iter()
            .any(|&n| mol.neighbors(n).filter(|&nb| is_hetero(nb)).count() >= 3)
        {
            fp.set(106);
        }

        // Key 107: [F,Cl,Br,I]~*(~*)~* = halogen on a trisubstituted (3+ bonds) atom
        if edges.iter().any(|&e| {
            let (u, v) = mol.edge_endpoints(e).unwrap();
            let check = |hal: petgraph::graph::NodeIndex, center: petgraph::graph::NodeIndex| {
                is_halogen(hal) && mol.neighbors(center).count() >= 3
            };
            check(u, v) || check(v, u)
        }) {
            fp.set(107);
        }

        // Key 108: [CH3]~*~*~*~[CH2]~* = 6-atom path; CH2 must have an outer heavy neighbor.
        // has_path5(&is_ch3, &is_ch2) finds 5-atom paths; the trailing * requires CH2 degree ≥ 2.
        if nodes.iter().filter(|&&n| is_ch3(n)).any(|&start| {
            // DFS 4 hops from CH3, check terminal is CH2 with outer neighbor
            let n = mol.node_count();
            if n <= 64 {
                let init_mask = 1u64 << start.index();
                let mut stack: Vec<(NodeIndex, usize, u64)> = vec![(start, 4, init_mask)];
                let mut found = false;
                while let Some((cur, rem, vis)) = stack.pop() {
                    if rem == 0 {
                        if is_ch2(cur) && mol.neighbors(cur).any(|nb| vis >> nb.index() & 1 == 0) {
                            found = true;
                            break;
                        }
                        continue;
                    }
                    for nb in mol.neighbors(cur) {
                        let bit = 1u64 << nb.index();
                        if vis & bit == 0 {
                            stack.push((nb, rem - 1, vis | bit));
                        }
                    }
                }
                found
            } else {
                has_simple_path(mol, start, 4, &|end| {
                    is_ch2(end) && mol.neighbors(end).count() >= 2
                })
            }
        }) {
            fp.set(108);
        }

        // Key 109: *~[CH2]~[#8] = CH2 bonded to O
        if has_path3(&is_heavy, &is_ch2, &is_oxygen) {
            fp.set(109);
        }

        // Key 110: [#7]~[#6]~[#8] = N-C-O
        if has_path3(&is_nitrogen, &is_carbon, &is_oxygen) {
            fp.set(110);
        }

        // Key 111: [#7]~*~[CH2]~* = N-*-CH2-*
        if has_path4(&is_nitrogen, &is_heavy, &is_ch2, &is_heavy) {
            fp.set(111);
        }

        // Key 112: *~*(~*)(~*)~* = atom with 4+ connections (quaternary center)
        if nodes.iter().any(|&n| mol.neighbors(n).count() >= 4) {
            fp.set(112);
        }

        // Key 113: [#8]!:*:* = O bonded (non-aromatically) to an aromatic ring
        if edges.iter().any(|&e| {
            if bond_is_aromatic(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            (is_oxygen(u) && is_aromatic_atom(v)) || (is_oxygen(v) && is_aromatic_atom(u))
        }) {
            fp.set(113);
        }

        // Key 114: [CH3]~[CH2]~* = CH3-CH2-*
        if has_path3(&is_ch3, &is_ch2, &is_heavy) {
            fp.set(114);
        }

        // Key 115: [CH3]~*~[CH2]~* = CH3-*-CH2-*
        if has_path4(&is_ch3, &is_heavy, &is_ch2, &is_heavy) {
            fp.set(115);
        }

        // Key 116: [$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]
        // Pattern 1: CH3~A~B~CH2~E (5 atoms; CH2 must have outer neighbor E)
        // Pattern 2: CH3 attached to a 3-membered ring containing CH2
        {
            let found116_p1 = edges.iter().any(|&e| {
                let (u, v) = mol.edge_endpoints(e).unwrap();
                let check = |b: NodeIndex, c: NodeIndex| -> bool {
                    if !is_heavy(b) || !is_heavy(c) {
                        return false;
                    }
                    for a in mol.neighbors(b) {
                        if a == c || !is_ch3(a) {
                            continue;
                        }
                        for d in mol.neighbors(c) {
                            if d == b || d == a || !is_ch2(d) {
                                continue;
                            }
                            // Trailing ~*: CH2(d) needs at least one heavy neighbor beyond c
                            if mol.neighbors(d).any(|ex| ex != c) {
                                return true;
                            }
                        }
                    }
                    false
                };
                check(u, v) || check(v, u)
            });
            // Pattern 2: [CH3]~*1~*~[CH2]1 = CH3 bonded to a 3-ring that contains CH2
            let found116_p2 = !found116_p1
                && nodes.iter().any(|&ch3n| {
                    if !is_ch3(ch3n) || !in_ring(ch3n) || !mol[ch3n].ring_sizes.contains(&3) {
                        return false;
                    }
                    mol.neighbors(ch3n).any(|ring_a| {
                        mol.neighbors(ring_a).any(|ring_b| {
                            ring_b != ch3n
                                && is_ch2(ring_b)
                                && mol.find_edge(ring_b, ch3n).is_some()
                        })
                    })
                });
            if found116_p1 || found116_p2 {
                fp.set(116);
            }
        }

        // Key 117: [#7]~*~[#8] = N-*-O (two bonds)
        if has_path3(&is_nitrogen, &is_heavy, &is_oxygen) {
            fp.set(117);
        }

        // Key 118: [$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)] count > 1
        // Single-atom recursive SMARTS: count atoms matching either alternative.
        // Atom `a` matches if:
        //   (1) a~CH2(b)~CH2(c)~heavy(d), all distinct (injective) — linear path
        //   (2) a is in a 3-membered ring with two CH2 neighbors: a-b(CH2)-c(CH2)-a
        // Key fires when count > 1.
        {
            let count_118 = nodes
                .iter()
                .filter(|&&a| {
                    // Condition 1: linear path a → CH2(b) → CH2(c, c≠a) → d (d≠b, d≠a)
                    let cond1 = mol.neighbors(a).any(|b| {
                        if !is_ch2(b) {
                            return false;
                        }
                        mol.neighbors(b).filter(|&c| c != a).any(|c| {
                            if !is_ch2(c) {
                                return false;
                            }
                            mol.neighbors(c).any(|d| d != b && d != a)
                        })
                    });
                    if cond1 {
                        return true;
                    }
                    // Condition 2: 3-ring a-b(CH2)-c(CH2)-a
                    mol.neighbors(a).any(|b| {
                        if !is_ch2(b) {
                            return false;
                        }
                        mol.neighbors(b)
                            .filter(|&c| c != a)
                            .any(|c| is_ch2(c) && mol.neighbors(c).any(|d| d == a))
                    })
                })
                .count();
            if count_118 > 1 {
                fp.set(118);
            }
        }

        // Key 119: [#7]=* = N in double bond
        if edges.iter().any(|&e| {
            bond_is_double(e) && {
                let (u, v) = mol.edge_endpoints(e).unwrap();
                is_nitrogen(u) || is_nitrogen(v)
            }
        }) {
            fp.set(119);
        }

        // Key 120: [!#6;R] min_count=1 → at least 2 non-C ring atoms
        if n_hetero_in_ring >= 2 {
            fp.set(120);
        }

        // Key 121: [#7;R] = N in ring
        if nodes.iter().any(|&n| is_nitrogen(n) && in_ring(n)) {
            fp.set(121);
        }

        // Key 122: *~[#7](~*)~* = N with 3+ bonds (tertiary N, any kind)
        if nodes
            .iter()
            .filter(|&&n| is_nitrogen(n))
            .any(|&n| mol.neighbors(n).count() >= 3)
        {
            fp.set(122);
        }

        // Key 123: [#8]~[#6]~[#8] = O-C-O
        if has_path3(&is_oxygen, &is_carbon, &is_oxygen) {
            fp.set(123);
        }

        // Key 124: [!#6;!#1]~[!#6;!#1] = heteroatom-heteroatom bond (any)
        if has_bond(&is_hetero, &is_hetero) {
            fp.set(124);
        }

        // Key 125: ('?', 0) = RDKit CalcNumAromaticRings > 1.
        // Use SSSR to avoid counting non-minimal "outer" rings in fused systems.
        // A ring is aromatic iff ALL consecutive atom-pair bonds in the ring are aromatic.
        {
            let sssr = molprint_core::ring::find_sssr(mol);
            let n_arom_rings = sssr
                .iter()
                .filter(|ring| {
                    let n = ring.len();
                    (0..n).all(|i| {
                        let a = ring[i];
                        let b = ring[(i + 1) % n];
                        mol.find_edge(a, b)
                            .map(|e| mol[e] == BondType::Aromatic)
                            .unwrap_or(false)
                    })
                })
                .count();
            if n_arom_rings > 1 {
                fp.set(125);
            }
        }

        // Key 126: *!@[#8]!@* = O with two non-ring bonds (ether O not in ring)
        if nodes
            .iter()
            .filter(|&&n| is_oxygen(n) && !in_ring(n))
            .any(|&o| mol.neighbors(o).count() >= 2)
        {
            fp.set(126);
        }

        // Key 127: *@*!@[#8] min_count=1 → count of (ring_neighbor, ring_atom, O) triples > 1
        // Each ring atom B bonded to O via non-ring bond contributes N_ring_bonds(B) triples.
        // Since ring atoms always have >=2 ring bonds, any ring atom with an O bond gives count>=2.
        {
            let count_127: usize = nodes
                .iter()
                .filter(|&&b| in_ring(b))
                .map(|&b| {
                    let n_ring_bonds = edges
                        .iter()
                        .filter(|&&e2| bond_in_ring(e2))
                        .filter(|&&e2| {
                            let (u, v) = mol.edge_endpoints(e2).unwrap();
                            u == b || v == b
                        })
                        .count();
                    let n_o_nonring = edges
                        .iter()
                        .filter(|&&e2| !bond_in_ring(e2))
                        .filter(|&&e2| {
                            let (u, v) = mol.edge_endpoints(e2).unwrap();
                            (u == b && is_oxygen(v)) || (v == b && is_oxygen(u))
                        })
                        .count();
                    n_ring_bonds * n_o_nonring
                })
                .sum();
            if count_127 > 1 {
                fp.set(127);
            }
        }

        // Key 128: [$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),
        //          $(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)] min_count=0
        // Pattern 1: acyclic *~CH2~A~B~C~CH2~* (4 hops, 7 distinct atoms)
        // Pattern 2: 6-ring where anchor's two ring-adjacent atoms are both ring-CH2
        // Pattern 3: CH2 attached to ring-atom A in 4-ring which has ring-CH2 ring-neighbor
        // Pattern 4: CH2, intermediate atom, ring-atom A in 3-ring with ring-CH2 ring-neighbor
        {
            let nc = mol.node_count();

            // Pattern 1: *~[CH2_a]~A~B~C~[CH2_b]~* (7 distinct atoms)
            let found128_p1 = if nc <= 64 {
                nodes.iter().filter(|&&n| is_ch2(n)).any(|&ch2_a| {
                    let init_mask = 1u64 << ch2_a.index();
                    let mut stack: Vec<(NodeIndex, usize, u64, NodeIndex)> = mol
                        .neighbors(ch2_a)
                        .map(|first| {
                            let m = init_mask | (1u64 << first.index());
                            (first, 3usize, m, first)
                        })
                        .collect();
                    while let Some((cur, rem, vis, first_hop)) = stack.pop() {
                        if rem == 0 {
                            if !is_ch2(cur) {
                                continue;
                            }
                            let mut found_pair = false;
                            'outer128: for a_out in mol.neighbors(ch2_a) {
                                if vis >> a_out.index() & 1 != 0 {
                                    continue;
                                }
                                for b_out in mol.neighbors(cur) {
                                    if vis >> b_out.index() & 1 != 0 {
                                        continue;
                                    }
                                    if a_out != b_out {
                                        found_pair = true;
                                        break 'outer128;
                                    }
                                }
                            }
                            if found_pair {
                                return true;
                            }
                            continue;
                        }
                        for nb in mol.neighbors(cur) {
                            let bit = 1u64 << nb.index();
                            if vis & bit == 0 {
                                stack.push((nb, rem - 1, vis | bit, first_hop));
                            }
                        }
                    }
                    false
                })
            } else {
                nodes.iter().filter(|&&n| is_ch2(n)).any(|&ch2_a| {
                    let mut stack: Vec<(
                        NodeIndex,
                        usize,
                        std::collections::HashSet<NodeIndex>,
                        NodeIndex,
                    )> = mol
                        .neighbors(ch2_a)
                        .map(|first| {
                            let mut vis = std::collections::HashSet::new();
                            vis.insert(ch2_a);
                            vis.insert(first);
                            (first, 3usize, vis, first)
                        })
                        .collect();
                    while let Some((cur, rem, vis, _first_hop)) = stack.pop() {
                        if rem == 0 {
                            if !is_ch2(cur) {
                                continue;
                            }
                            let mut found_pair = false;
                            'outer128b: for a_out in mol.neighbors(ch2_a) {
                                if vis.contains(&a_out) {
                                    continue;
                                }
                                for b_out in mol.neighbors(cur) {
                                    if vis.contains(&b_out) {
                                        continue;
                                    }
                                    if a_out != b_out {
                                        found_pair = true;
                                        break 'outer128b;
                                    }
                                }
                            }
                            if found_pair {
                                return true;
                            }
                            continue;
                        }
                        for nb in mol.neighbors(cur) {
                            if !vis.contains(&nb) {
                                let mut nv = vis.clone();
                                nv.insert(nb);
                                stack.push((nb, rem - 1, nv, ch2_a));
                            }
                        }
                    }
                    false
                })
            };

            // Pattern 2: 6-membered ring where anchor's two ring-adjacent atoms are both CH2
            // (ring-bond DFS; all traversed atoms are in rings, so is_ch2 suffices)
            let found128_p2 = !found128_p1 && {
                let mut found = false;
                'p2_128: for &start in &nodes {
                    if !in_ring(start) {
                        continue;
                    }
                    if nc <= 64 {
                        let init_vis: u64 = 1u64 << start.index();
                        let mut stack: Vec<(NodeIndex, usize, u64, NodeIndex)> = Vec::new();
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            stack.push((nb, 1usize, init_vis | (1u64 << nb.index()), nb));
                        }
                        while let Some((cur, depth, vis, first_step)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 5 && is_ch2(first_step) && is_ch2(cur) {
                                        found = true;
                                        break 'p2_128;
                                    }
                                    continue;
                                }
                                let bit = 1u64 << nb.index();
                                if vis & bit != 0 {
                                    continue;
                                }
                                if depth < 5 {
                                    stack.push((nb, depth + 1, vis | bit, first_step));
                                }
                            }
                            if found {
                                break;
                            }
                        }
                    } else {
                        let mut init_vis = std::collections::HashSet::new();
                        init_vis.insert(start);
                        let mut stack: Vec<(
                            NodeIndex,
                            usize,
                            std::collections::HashSet<NodeIndex>,
                            NodeIndex,
                        )> = Vec::new();
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            let mut vis = init_vis.clone();
                            vis.insert(nb);
                            stack.push((nb, 1usize, vis, nb));
                        }
                        while let Some((cur, depth, vis, first_step)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 5 && is_ch2(first_step) && is_ch2(cur) {
                                        found = true;
                                        break 'p2_128;
                                    }
                                    continue;
                                }
                                if vis.contains(&nb) {
                                    continue;
                                }
                                if depth < 5 {
                                    let mut nv = vis.clone();
                                    nv.insert(nb);
                                    stack.push((nb, depth + 1, nv, first_step));
                                }
                            }
                            if found {
                                break;
                            }
                        }
                    }
                }
                found
            };

            // Pattern 3: *~[CH2]~A(in 4-ring)~...~D(ring-CH2, in same 4-ring)~back to A
            // Traverse 4-ring A-B-C-D-A; find specific members {B,C,D}; E must not be in {B,C,D}.
            // E can be ring-bonded to A from a different ring (e.g. spiro atom).
            let found128_p3 = !found128_p1
                && !found128_p2
                && nodes.iter().any(|&a| {
                    if !in_ring(a) || !mol[a].ring_sizes.contains(&4) {
                        return false;
                    }
                    let ring_nbs_a: Vec<NodeIndex> = mol
                        .edges(a)
                        .filter(|er| ring_bond_set.contains(&er.id()))
                        .map(|er| {
                            if er.source() == a {
                                er.target()
                            } else {
                                er.source()
                            }
                        })
                        .collect();
                    // Find a specific 4-ring A-B-C-D-A where D is ring-CH2; collect {B,C,D}.
                    let ring4 = ring_nbs_a.iter().find_map(|&b| {
                        mol.edges(b)
                            .filter(|er| ring_bond_set.contains(&er.id()))
                            .map(|er| {
                                if er.source() == b {
                                    er.target()
                                } else {
                                    er.source()
                                }
                            })
                            .filter(|&c| c != a)
                            .find_map(|c| {
                                mol.edges(c)
                                    .filter(|er| ring_bond_set.contains(&er.id()))
                                    .map(|er| {
                                        if er.source() == c {
                                            er.target()
                                        } else {
                                            er.source()
                                        }
                                    })
                                    .find(|&d| {
                                        d != b && ring_nbs_a.contains(&d) && is_ch2(d) && in_ring(d)
                                    })
                                    .map(|d| (b, c, d))
                            })
                    });
                    let (rb, rc, rd) = match ring4 {
                        None => return false,
                        Some(x) => x,
                    };
                    // E is any CH2 neighbor of A NOT in the 4-ring {B,C,D}, with outer neighbor
                    mol.neighbors(a).any(|e| {
                        e != rb
                            && e != rc
                            && e != rd
                            && is_ch2(e)
                            && mol.neighbors(e).any(|x| x != a)
                    })
                });

            // Pattern 4: *~[CH2]~X~A(in 3-ring)~...~C(ring-CH2, in same 3-ring)~back to A
            // 3-ring: triangle A-B-C-A; find {B,C}; X must not be B or C.
            // X can be ring-bonded to A from a different ring.
            let found128_p4 = !found128_p1
                && !found128_p2
                && !found128_p3
                && nodes.iter().any(|&a| {
                    if !in_ring(a) || !mol[a].ring_sizes.contains(&3) {
                        return false;
                    }
                    let ring_nbs_a: Vec<NodeIndex> = mol
                        .edges(a)
                        .filter(|er| ring_bond_set.contains(&er.id()))
                        .map(|er| {
                            if er.source() == a {
                                er.target()
                            } else {
                                er.source()
                            }
                        })
                        .collect();
                    // Find a specific 3-ring A-B-C-A where C is ring-CH2; collect (B, C).
                    let ring3 = ring_nbs_a.iter().find_map(|&b| {
                        mol.edges(b)
                            .filter(|er| ring_bond_set.contains(&er.id()))
                            .map(|er| {
                                if er.source() == b {
                                    er.target()
                                } else {
                                    er.source()
                                }
                            })
                            .find(|&c| c != a && ring_nbs_a.contains(&c) && is_ch2(c) && in_ring(c))
                            .map(|c| (b, c))
                    });
                    let (rb, rc) = match ring3 {
                        None => return false,
                        Some(x) => x,
                    };
                    // X is any neighbor of A NOT in the 3-ring {B,C}; X has CH2 neighbor E (≠A)
                    mol.neighbors(a).any(|x| {
                        x != rb
                            && x != rc
                            && mol.neighbors(x).any(|e| {
                                e != a && is_ch2(e) && mol.neighbors(e).any(|outer| outer != x)
                            })
                    })
                });

            if found128_p1 || found128_p2 || found128_p3 || found128_p4 {
                fp.set(128);
            }
        }

        // Key 129: [$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),
        //          $(*~[CH2]~[R]1@[R]@[CH2;R]1)] min_count=0
        // Pattern 1: *~CH2~A~B~CH2~* (3 hops, 6 distinct atoms)
        // Pattern 2: 5-ring where anchor's two ring-adjacent atoms are both CH2
        // Pattern 3: CH2 attached to ring-atom A in 3-ring which has ring-CH2 ring-neighbor
        {
            let nc = mol.node_count();

            // Pattern 1: *~[CH2_a]~A~B~[CH2_b]~* (6 distinct atoms)
            let found129_p1 = edges.iter().any(|&bc| {
                let (b, c) = mol.edge_endpoints(bc).unwrap();
                let try_bc = |ch2_b: NodeIndex, mid_c: NodeIndex| -> bool {
                    if !is_ch2(ch2_b) {
                        return false;
                    }
                    for a in mol.neighbors(ch2_b) {
                        if a == mid_c {
                            continue;
                        }
                        for mid_d in mol.neighbors(mid_c) {
                            if mid_d == ch2_b || mid_d == a {
                                continue;
                            }
                            for e in mol.neighbors(mid_d) {
                                if e == mid_c || e == ch2_b || e == a {
                                    continue;
                                }
                                if !is_ch2(e) {
                                    continue;
                                }
                                if mol
                                    .neighbors(e)
                                    .any(|f| f != mid_d && f != mid_c && f != ch2_b && f != a)
                                {
                                    return true;
                                }
                            }
                        }
                    }
                    false
                };
                try_bc(b, c) || try_bc(c, b)
            });

            // Pattern 2: 5-membered ring where anchor's two ring-adjacent atoms are both CH2
            let found129_p2 = !found129_p1 && {
                let mut found = false;
                'p2_129: for &start in &nodes {
                    if !in_ring(start) {
                        continue;
                    }
                    if nc <= 64 {
                        let init_vis: u64 = 1u64 << start.index();
                        let mut stack: Vec<(NodeIndex, usize, u64, NodeIndex)> = Vec::new();
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            stack.push((nb, 1usize, init_vis | (1u64 << nb.index()), nb));
                        }
                        while let Some((cur, depth, vis, first_step)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 4 && is_ch2(first_step) && is_ch2(cur) {
                                        found = true;
                                        break 'p2_129;
                                    }
                                    continue;
                                }
                                let bit = 1u64 << nb.index();
                                if vis & bit != 0 {
                                    continue;
                                }
                                if depth < 4 {
                                    stack.push((nb, depth + 1, vis | bit, first_step));
                                }
                            }
                            if found {
                                break;
                            }
                        }
                    } else {
                        let mut init_vis = std::collections::HashSet::new();
                        init_vis.insert(start);
                        let mut stack: Vec<(
                            NodeIndex,
                            usize,
                            std::collections::HashSet<NodeIndex>,
                            NodeIndex,
                        )> = Vec::new();
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            let mut vis = init_vis.clone();
                            vis.insert(nb);
                            stack.push((nb, 1usize, vis, nb));
                        }
                        while let Some((cur, depth, vis, first_step)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 4 && is_ch2(first_step) && is_ch2(cur) {
                                        found = true;
                                        break 'p2_129;
                                    }
                                    continue;
                                }
                                if vis.contains(&nb) {
                                    continue;
                                }
                                if depth < 4 {
                                    let mut nv = vis.clone();
                                    nv.insert(nb);
                                    stack.push((nb, depth + 1, nv, first_step));
                                }
                            }
                            if found {
                                break;
                            }
                        }
                    }
                }
                found
            };

            // Pattern 3: *~[CH2]~A(in 3-ring) where A has ring-CH2 ring-neighbor.
            // Find {B,C} in the 3-ring; E must not be B or C (can be ring-bonded via diff ring).
            let found129_p3 = !found129_p1
                && !found129_p2
                && nodes.iter().any(|&a| {
                    if !in_ring(a) || !mol[a].ring_sizes.contains(&3) {
                        return false;
                    }
                    let ring_nbs_a: Vec<NodeIndex> = mol
                        .edges(a)
                        .filter(|er| ring_bond_set.contains(&er.id()))
                        .map(|er| {
                            if er.source() == a {
                                er.target()
                            } else {
                                er.source()
                            }
                        })
                        .collect();
                    // Find a specific 3-ring A-B-C-A where C is ring-CH2; collect (B, C).
                    let ring3 = ring_nbs_a.iter().find_map(|&b| {
                        mol.edges(b)
                            .filter(|er| ring_bond_set.contains(&er.id()))
                            .map(|er| {
                                if er.source() == b {
                                    er.target()
                                } else {
                                    er.source()
                                }
                            })
                            .find(|&c| c != a && ring_nbs_a.contains(&c) && is_ch2(c) && in_ring(c))
                            .map(|c| (b, c))
                    });
                    let (rb, rc) = match ring3 {
                        None => return false,
                        Some(x) => x,
                    };
                    // E is any CH2 neighbor of A NOT in the 3-ring {B,C}, with outer neighbor
                    mol.neighbors(a).any(|e| {
                        e != rb && e != rc && is_ch2(e) && mol.neighbors(e).any(|x| x != a)
                    })
                });

            if found129_p1 || found129_p2 || found129_p3 {
                fp.set(129);
            }
        }

        // Key 130: [!#6;!#1]~[!#6;!#1] min_count=1 → at least 2 heteroatom-heteroatom bonds
        {
            let count = count_bonded(&is_hetero, &is_hetero);
            if count >= 2 {
                fp.set(130);
            }
        }

        // Key 131: [!#6;!#1;!H0] min_count=1 → at least 2 H-bearing heteroatoms
        if n_hetero_with_h >= 2 {
            fp.set(131);
        }

        // Key 132: [#8]~*~[CH2]~* = O-*-CH2-*
        if has_path4(&is_oxygen, &is_heavy, &is_ch2, &is_heavy) {
            fp.set(132);
        }

        // Key 133: *@*!@[#7] = ring atom bonded to N via non-ring bond
        if edges.iter().any(|&e| {
            if bond_in_ring(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            (in_ring(u) && is_nitrogen(v)) || (in_ring(v) && is_nitrogen(u))
        }) {
            fp.set(133);
        }

        // Key 134: [F,Cl,Br,I] = any halogen
        if nodes.iter().any(|&n| is_halogen(n)) {
            fp.set(134);
        }

        // Key 135: [#7]!:*:* = N bonded (non-aromatically) to aromatic ring
        if edges.iter().any(|&e| {
            if bond_is_aromatic(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            (is_nitrogen(u) && is_aromatic_atom(v)) || (is_nitrogen(v) && is_aromatic_atom(u))
        }) {
            fp.set(135);
        }

        // Key 136: [#8]=* min_count=1 → at least 2 C=O (or any X=O)
        {
            let count = edges
                .iter()
                .filter(|&&e| {
                    bond_is_double(e) && {
                        let (u, v) = mol.edge_endpoints(e).unwrap();
                        is_oxygen(u) || is_oxygen(v)
                    }
                })
                .count();
            if count >= 2 {
                fp.set(136);
            }
        }

        // Key 137: [!C;!c;R] = non-carbon atom in ring (heteroatom in ring)
        if n_hetero_in_ring >= 1 {
            fp.set(137);
        }

        // Key 138: [!#6;!#1]~[CH2]~* min_count=1 → count of unique SMARTS matches (uniquify=True) > 1.
        // RDKit uniquify=True: two matches are the same iff they use the same unordered atom set.
        // So count = number of distinct 3-atom sets {hetero, ch2, trailing*}.
        {
            let mut unique138: std::collections::HashSet<(usize, usize, usize)> =
                std::collections::HashSet::new();
            for &c2 in &nodes {
                if !is_ch2(c2) {
                    continue;
                }
                let c2_nb: Vec<NodeIndex> = mol.neighbors(c2).collect();
                for &a1 in &c2_nb {
                    if !is_hetero(a1) {
                        continue;
                    }
                    for &a3 in &c2_nb {
                        if a3 == a1 {
                            continue;
                        }
                        let mut arr = [a1.index(), c2.index(), a3.index()];
                        arr.sort_unstable();
                        unique138.insert((arr[0], arr[1], arr[2]));
                    }
                }
            }
            if unique138.len() > 1 {
                fp.set(138);
            }
        }

        // Key 139: [O;!H0] = OH (O with H)
        if nodes.iter().any(|&n| is_oh(n)) {
            fp.set(139);
        }

        // Key 140: [#8] min_count=3 → at least 4 O atoms
        if n_oxygen >= 4 {
            fp.set(140);
        }

        // Key 141: [CH3] min_count=2 → at least 3 CH3 groups
        if n_ch3 >= 3 {
            fp.set(141);
        }

        // Key 142: [#7] min_count=1 → at least 2 N atoms
        if n_nitrogen >= 2 {
            fp.set(142);
        }

        // Key 143: *@*!@[#8] = ring atom bonded to O via non-ring bond (same as 127 but min_count=0)
        if edges.iter().any(|&e| {
            if bond_in_ring(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            (in_ring(u) && is_oxygen(v)) || (in_ring(v) && is_oxygen(u))
        }) {
            fp.set(143);
        }

        // Key 144: *!:*:*!:* = aromatic bond where BOTH endpoints have a non-aromatic bond
        // i.e., path: (non-arom bond)-[arom atom]:(arom atom)-(non-arom bond)
        if edges.iter().any(|&e| {
            if !bond_is_aromatic(e) {
                return false;
            }
            let (u, v) = mol.edge_endpoints(e).unwrap();
            let has_non_arom =
                |n: petgraph::graph::NodeIndex| mol.edges(n).any(|er| !bond_is_aromatic(er.id()));
            has_non_arom(u) && has_non_arom(v)
        }) {
            fp.set(144);
        }

        // Key 145: *1~*~*~*~*~*~1 min_count=1 → at least 2 six-membered rings (incl. non-SSSR).
        // Count distinct 6-cycles in the ring-bond subgraph. Each ring counted exactly once by
        // requiring start = min-index node AND first_step.index() < cur.index() (canonical direction).
        {
            let n = mol.node_count();
            let mut count_6 = 0usize;
            for &start in &nodes {
                if !in_ring(start) {
                    continue;
                }
                if n <= 64 {
                    let init_vis: u64 = 1u64 << start.index();
                    // Stack: (cur, depth, vis, first_step)
                    let mut stack: Vec<(NodeIndex, usize, u64, NodeIndex)> = Vec::new();
                    for er in mol.edges(start) {
                        if !ring_bond_set.contains(&er.id()) {
                            continue;
                        }
                        let nb = if er.source() == start {
                            er.target()
                        } else {
                            er.source()
                        };
                        if nb.index() > start.index() {
                            stack.push((nb, 1, init_vis | (1u64 << nb.index()), nb));
                        }
                    }
                    while let Some((cur, depth, vis, first_step)) = stack.pop() {
                        for er in mol.edges(cur) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == cur {
                                er.target()
                            } else {
                                er.source()
                            };
                            if nb == start {
                                // Canonical direction: first_step < last-node-before-start
                                if depth == 5 && first_step.index() < cur.index() {
                                    count_6 += 1;
                                }
                                continue;
                            }
                            if nb.index() <= start.index() {
                                continue;
                            }
                            let bit = 1u64 << nb.index();
                            if vis & bit != 0 {
                                continue;
                            }
                            if depth < 5 {
                                stack.push((nb, depth + 1, vis | bit, first_step));
                            }
                        }
                    }
                } else {
                    let mut stack: Vec<(
                        NodeIndex,
                        usize,
                        std::collections::HashSet<NodeIndex>,
                        NodeIndex,
                    )> = Vec::new();
                    let mut init_vis = std::collections::HashSet::new();
                    init_vis.insert(start);
                    for er in mol.edges(start) {
                        if !ring_bond_set.contains(&er.id()) {
                            continue;
                        }
                        let nb = if er.source() == start {
                            er.target()
                        } else {
                            er.source()
                        };
                        if nb.index() > start.index() {
                            let mut vis = init_vis.clone();
                            vis.insert(nb);
                            stack.push((nb, 1, vis, nb));
                        }
                    }
                    while let Some((cur, depth, vis, first_step)) = stack.pop() {
                        for er in mol.edges(cur) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == cur {
                                er.target()
                            } else {
                                er.source()
                            };
                            if nb == start {
                                if depth == 5 && first_step.index() < cur.index() {
                                    count_6 += 1;
                                }
                                continue;
                            }
                            if nb.index() <= start.index() {
                                continue;
                            }
                            if vis.contains(&nb) {
                                continue;
                            }
                            if depth < 5 {
                                let mut nv = vis.clone();
                                nv.insert(nb);
                                stack.push((nb, depth + 1, nv, first_step));
                            }
                        }
                    }
                }
                if count_6 >= 2 {
                    break;
                }
            }
            if count_6 >= 2 {
                fp.set(145);
            }
        }

        // Key 146: [#8] min_count=2 → at least 3 O atoms
        if n_oxygen >= 3 {
            fp.set(146);
        }

        // Key 147: [$(*~[CH2]~[CH2]~*),...] = any CH2-CH2 sequence
        if edges.iter().any(|&e| {
            let (u, v) = mol.edge_endpoints(e).unwrap();
            is_ch2(u) && is_ch2(v)
        }) {
            fp.set(147);
        }

        // Key 148: *~[!#6;!#1](~*)~* = heteroatom with 3+ bonds
        if nodes
            .iter()
            .filter(|&&n| is_hetero(n))
            .any(|&n| mol.neighbors(n).count() >= 3)
        {
            fp.set(148);
        }

        // Key 149: [C;H3,H4] min_count=1 → at least 2 terminal C
        if n_ch3_or_ch4 >= 2 {
            fp.set(149);
        }

        // Key 150: *!@*@*!@* = ring bond (B@C) where both B and C have a non-ring bond
        if edges.iter().any(|&e| {
            if !bond_in_ring(e) {
                return false;
            }
            let (b, c) = mol.edge_endpoints(e).unwrap();
            let has_non_ring_bond = |n: petgraph::graph::NodeIndex| {
                edges.iter().any(|&e2| {
                    if bond_in_ring(e2) {
                        return false;
                    }
                    let (u, v) = mol.edge_endpoints(e2).unwrap();
                    u == n || v == n
                })
            };
            has_non_ring_bond(b) && has_non_ring_bond(c)
        }) {
            fp.set(150);
        }

        // Key 151: [#7;!H0] = N with at least one H
        if nodes.iter().any(|&n| is_n_with_h(n)) {
            fp.set(151);
        }

        // Key 152: [#8]~[#6](~[#6])~[#6] = O-C(-C)-C (O bonded to C that has 2+ C neighbors)
        if nodes.iter().filter(|&&n| is_carbon(n)).any(|&c| {
            let o_count = mol.neighbors(c).filter(|&nb| is_oxygen(nb)).count();
            let c_count = mol.neighbors(c).filter(|&nb| is_carbon(nb)).count();
            o_count >= 1 && c_count >= 2
        }) {
            fp.set(152);
        }

        // Key 153: [!#6;!#1]~[CH2]~* = heteroatom-CH2-*
        if has_path3(&is_hetero, &is_ch2, &is_heavy) {
            fp.set(153);
        }

        // Key 154: [#6]=[#8] = C=O
        if has_double_bond(&is_carbon, &is_oxygen) {
            fp.set(154);
        }

        // Key 155: *!@[CH2]!@* = non-ring CH2 with TWO non-ring bonds to heavy atoms
        // CH2 must have ≥2 heavy-atom neighbors (degree ≥ 2), all via non-ring bonds.
        // If CH2 is not in a ring, all its bonds are non-ring bonds automatically.
        if nodes
            .iter()
            .any(|&n| is_ch2(n) && !in_ring(n) && mol.neighbors(n).count() >= 2)
        {
            fp.set(155);
        }

        // Key 156: [#7]~*(~*)~* = N bonded to an atom with heavy-atom degree >= 3.
        // Iterate over N nodes to avoid edge-direction bug when both endpoints are N.
        if nodes
            .iter()
            .any(|&n| is_nitrogen(n) && mol.neighbors(n).any(|nb| mol.neighbors(nb).count() >= 3))
        {
            fp.set(156);
        }

        // Key 157: [#6]-[#8] = C-O single bond
        if has_single_bond(&is_carbon, &is_oxygen) {
            fp.set(157);
        }

        // Key 158: [#6]-[#7] = C-N single bond
        if has_single_bond(&is_carbon, &is_nitrogen) {
            fp.set(158);
        }

        // Key 159: [#8] min_count=1 → at least 2 O atoms
        if n_oxygen >= 2 {
            fp.set(159);
        }

        // Key 160: [C;H3,H4] = at least 1 terminal C (methyl or methane)
        if n_ch3_or_ch4 >= 1 {
            fp.set(160);
        }

        // Key 161: [#7] = at least 1 N
        if n_nitrogen >= 1 {
            fp.set(161);
        }

        // Key 162: a = any aromatic atom
        if nodes.iter().any(|&n| is_aromatic_atom(n)) {
            fp.set(162);
        }

        // Key 163: *1~*~*~*~*~*~1 = any 6-membered ring (incl. non-SSSR outer rings)
        {
            let has_6_ring = has_ring_size(6) || {
                let n = mol.node_count();
                let mut found = false;
                'dfs163: for &start in &nodes {
                    if !in_ring(start) {
                        continue;
                    }
                    if n <= 64 {
                        let init_vis: u64 = 1u64 << start.index();
                        let mut stack: Vec<(NodeIndex, usize, u64)> = Vec::new();
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            stack.push((nb, 1, init_vis | (1u64 << nb.index())));
                        }
                        while let Some((cur, depth, vis)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 5 {
                                        found = true;
                                        break 'dfs163;
                                    }
                                    continue;
                                }
                                let bit = 1u64 << nb.index();
                                if vis & bit != 0 {
                                    continue;
                                }
                                if depth < 5 {
                                    stack.push((nb, depth + 1, vis | bit));
                                }
                            }
                        }
                    } else {
                        let mut stack: Vec<(
                            NodeIndex,
                            usize,
                            std::collections::HashSet<NodeIndex>,
                        )> = Vec::new();
                        let mut init_vis = std::collections::HashSet::new();
                        init_vis.insert(start);
                        for er in mol.edges(start) {
                            if !ring_bond_set.contains(&er.id()) {
                                continue;
                            }
                            let nb = if er.source() == start {
                                er.target()
                            } else {
                                er.source()
                            };
                            let mut vis = init_vis.clone();
                            vis.insert(nb);
                            stack.push((nb, 1, vis));
                        }
                        while let Some((cur, depth, vis)) = stack.pop() {
                            for er in mol.edges(cur) {
                                if !ring_bond_set.contains(&er.id()) {
                                    continue;
                                }
                                let nb = if er.source() == cur {
                                    er.target()
                                } else {
                                    er.source()
                                };
                                if nb == start {
                                    if depth == 5 {
                                        found = true;
                                        break 'dfs163;
                                    }
                                    continue;
                                }
                                if vis.contains(&nb) {
                                    continue;
                                }
                                if depth < 5 {
                                    let mut nv = vis.clone();
                                    nv.insert(nb);
                                    stack.push((nb, depth + 1, nv));
                                }
                            }
                        }
                    }
                }
                found
            };
            if has_6_ring {
                fp.set(163);
            }
        }

        // Key 164: [#8] = at least 1 O
        if n_oxygen >= 1 {
            fp.set(164);
        }

        // Key 165: [R] = any ring atom
        if n_in_ring >= 1 {
            fp.set(165);
        }

        // Key 166: `?` always 0

        fp
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molprint_core::smiles::parse_smiles;

    fn parse(s: &str) -> MolGraph {
        parse_smiles(s).expect(s)
    }

    fn bits(s: &str) -> Vec<usize> {
        let fp = Maccs166::new().fingerprint(&parse(s));
        let mut v: Vec<usize> = (0..166).filter(|&i| fp.get(i)).collect();
        v.sort_unstable();
        v
    }

    #[test]
    fn benzene_bits() {
        let b = bits("c1ccccc1");
        // Should have: 163 (6-ring), 165 (ring atom), 162 (aromatic)
        assert!(b.contains(&163), "missing 6-ring: {b:?}");
        assert!(b.contains(&165), "missing ring: {b:?}");
        assert!(b.contains(&162), "missing aromatic: {b:?}");
    }

    #[test]
    fn ethanol_bits() {
        let b = bits("CCO");
        // O present (164), C-O (157), OH (139)
        assert!(b.contains(&164), "missing O: {b:?}");
        assert!(b.contains(&157), "missing C-O: {b:?}");
        assert!(b.contains(&139), "missing OH: {b:?}");
        // No ring
        assert!(!b.contains(&165), "spurious ring: {b:?}");
    }

    #[test]
    fn methane_bits() {
        let b = bits("C");
        // [C;H3,H4] (160), no N, no O, no ring
        assert!(b.contains(&160), "missing CH4 (160): {b:?}");
        assert!(!b.contains(&161), "spurious N: {b:?}");
        assert!(!b.contains(&164), "spurious O: {b:?}");
        assert!(!b.contains(&165), "spurious ring: {b:?}");
    }

    #[test]
    fn pyridine_bits() {
        let b = bits("c1ccncc1");
        // aromatic N-ring
        assert!(b.contains(&162), "missing aromatic: {b:?}");
        assert!(b.contains(&163), "missing 6-ring: {b:?}");
        assert!(b.contains(&161), "missing N: {b:?}");
        assert!(b.contains(&121), "missing N in ring: {b:?}");
        assert!(b.contains(&98), "missing hetero 6-ring: {b:?}");
        assert!(b.contains(&65), "missing c:n: {b:?}");
    }

    #[test]
    fn nitromethane_bits() {
        // C[N+](=O)[O-] — charged, N=O, NO2
        let b = bits("C[N+](=O)[O-]");
        assert!(b.contains(&49), "missing charged: {b:?}");
        assert!(b.contains(&56), "missing nitro O-N(-O)-C: {b:?}");
    }
}
