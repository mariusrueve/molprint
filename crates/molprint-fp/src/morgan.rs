use crate::bitvec::FingerprintBits;
use crate::traits::Fingerprinter;
use molprint_core::mol::graph::MolGraphExt;
use molprint_core::MolGraph;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef as _;

/// Deterministic hash combine (boost::hash_combine style).
/// AHasher::default() uses a randomized seed and must NOT be used here.
fn hash_combine(h: u64, val: u64) -> u64 {
    h ^ val
        .wrapping_add(0x9e3779b97f4a7c15)
        .wrapping_add(h << 6)
        .wrapping_add(h >> 2)
}

fn hash_u64(val: u64) -> u64 {
    // FNV-1a inspired 64-bit hash
    let mut h: u64 = 0xcbf29ce484222325;
    for byte in val.to_le_bytes() {
        h ^= byte as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn combine_hashes(base: u64, vals: &[(u8, u32)]) -> u32 {
    let mut h = hash_u64(base);
    for &(bond, nb) in vals {
        h = hash_combine(h, bond as u64);
        h = hash_combine(h, nb as u64);
    }
    h as u32
}

pub struct Morgan {
    pub radius: u32,
    pub nbits: usize,
    pub use_bond_types: bool,
}

impl Morgan {
    pub fn new(radius: u32, nbits: usize) -> Self {
        Self {
            radius,
            nbits,
            use_bond_types: true,
        }
    }

    /// Compute atom invariant hash (radius-0 Morgan). Must be deterministic across runs.
    fn atom_invariant(&self, mol: &MolGraph, idx: NodeIndex) -> u32 {
        let atom = mol.atom(idx);
        let packed = ((atom.element.atomic_number() as u64) << 48)
            | ((mol.degree(idx) as u64) << 32)
            | ((atom.h_count as u64) << 16)
            | (((atom.charge as i64 + 128) as u64) << 8)
            | (atom.in_ring as u64);
        hash_u64(packed) as u32
    }

    /// One Morgan iteration: update hashes based on neighborhood.
    fn iterate(&self, mol: &MolGraph, prev_hashes: &[u32], _radius: u32) -> Vec<u32> {
        let n = mol.node_count();
        let mut new_hashes = vec![0u32; n];

        for i in 0..n {
            let idx = NodeIndex::new(i);
            let mut neighbors: Vec<(u8, u32)> = mol
                .edges(idx)
                .map(|e| {
                    let nb = if e.source() == idx {
                        e.target()
                    } else {
                        e.source()
                    };
                    let bond_order = if self.use_bond_types {
                        e.weight().order_x10()
                    } else {
                        10 // treat all as single
                    };
                    (bond_order, prev_hashes[nb.index()])
                })
                .collect();

            // Sort for deterministic ordering
            neighbors.sort_unstable();

            // If no neighbors, the environment doesn't change — keep same hash.
            if neighbors.is_empty() {
                new_hashes[i] = prev_hashes[i];
            } else {
                new_hashes[i] = combine_hashes(prev_hashes[i] as u64, &neighbors);
            }
        }

        new_hashes
    }
}

impl Fingerprinter for Morgan {
    fn fingerprint(&self, mol: &MolGraph) -> FingerprintBits {
        let n = mol.node_count();
        let mut fp = FingerprintBits::new(self.nbits);

        if n == 0 {
            return fp;
        }

        let mut current_hashes: Vec<u32> = (0..n)
            .map(|i| self.atom_invariant(mol, NodeIndex::new(i)))
            .collect();

        let mut all_hashes: std::collections::HashSet<u32> = std::collections::HashSet::new();

        for &h in &current_hashes {
            all_hashes.insert(h);
        }

        for r in 1..=self.radius {
            current_hashes = self.iterate(mol, &current_hashes, r);
            for &h in &current_hashes {
                all_hashes.insert(h);
            }
        }

        for h in all_hashes {
            fp.set(h as usize % self.nbits);
        }

        fp
    }

    fn name(&self) -> &str {
        "Morgan"
    }

    fn nbits(&self) -> usize {
        self.nbits
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molprint_core::smiles::parse_smiles;

    #[test]
    fn test_methane_fingerprint() {
        let mol = parse_smiles("C").unwrap();
        let morgan = Morgan::new(2, 2048);
        let fp = morgan.fingerprint(&mol);
        assert!(fp.count_ones() > 0);
        assert!(fp.count_ones() < 10);
    }

    #[test]
    fn test_same_molecule_same_fp() {
        let mol1 = parse_smiles("CCO").unwrap();
        let mol2 = parse_smiles("OCC").unwrap();
        let morgan = Morgan::new(2, 2048);
        let fp1 = morgan.fingerprint(&mol1);
        let fp2 = morgan.fingerprint(&mol2);
        assert_eq!(
            fp1, fp2,
            "same molecule must produce identical fingerprints"
        );
    }

    #[test]
    fn test_different_molecules_different_fp() {
        let mol1 = parse_smiles("CCO").unwrap();
        let mol2 = parse_smiles("CCCO").unwrap();
        let morgan = Morgan::new(2, 2048);
        let fp1 = morgan.fingerprint(&mol1);
        let fp2 = morgan.fingerprint(&mol2);
        assert_ne!(fp1, fp2);
    }

    #[test]
    fn test_benzene_more_bits_than_methane() {
        let methane = parse_smiles("C").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let morgan = Morgan::new(2, 2048);
        let fp_m = morgan.fingerprint(&methane);
        let fp_b = morgan.fingerprint(&benzene);
        assert!(fp_b.count_ones() > fp_m.count_ones());
    }

    #[test]
    fn test_radius_affects_fp() {
        let mol = parse_smiles("c1ccccc1CC").unwrap();
        let m1 = Morgan::new(1, 2048);
        let m2 = Morgan::new(3, 2048);
        let fp1 = m1.fingerprint(&mol);
        let fp2 = m2.fingerprint(&mol);
        assert!(fp2.count_ones() >= fp1.count_ones());
    }
}
