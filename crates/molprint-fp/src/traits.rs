use crate::bitvec::FingerprintBits;
use molprint_core::MolGraph;

/// Trait for all fingerprint algorithms.
pub trait Fingerprinter: Send + Sync {
    /// Compute fingerprint for a molecule.
    fn fingerprint(&self, mol: &MolGraph) -> FingerprintBits;

    /// Name of this fingerprint type (e.g., "Morgan", "MACCS166").
    fn name(&self) -> &str;

    /// Number of bits in the output fingerprint.
    fn nbits(&self) -> usize;
}
