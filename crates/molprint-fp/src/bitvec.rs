/// Cache-line-aligned bitvector for molecular fingerprints.
#[derive(Debug, Clone, PartialEq)]
pub struct FingerprintBits {
    words: Vec<u64>,
    nbits: usize,
}

impl FingerprintBits {
    /// Create a zero-initialized fingerprint of `nbits` bits.
    pub fn new(nbits: usize) -> Self {
        let nwords = nbits.div_ceil(64);
        Self {
            words: vec![0u64; nwords],
            nbits,
        }
    }

    /// Set bit at `bit_index` (already folded to [0, nbits)).
    pub fn set(&mut self, bit_index: usize) {
        let bit = bit_index % self.nbits;
        self.words[bit / 64] |= 1u64 << (bit % 64);
    }

    /// Get bit at `bit_index`.
    pub fn get(&self, bit_index: usize) -> bool {
        let bit = bit_index % self.nbits;
        (self.words[bit / 64] >> (bit % 64)) & 1 == 1
    }

    /// Count set bits (popcount).
    pub fn count_ones(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }

    /// Bitwise AND.
    pub fn and(&self, other: &FingerprintBits) -> FingerprintBits {
        let words = self
            .words
            .iter()
            .zip(other.words.iter())
            .map(|(a, b)| a & b)
            .collect();
        FingerprintBits {
            words,
            nbits: self.nbits,
        }
    }

    /// Bitwise OR.
    pub fn or(&self, other: &FingerprintBits) -> FingerprintBits {
        let words = self
            .words
            .iter()
            .zip(other.words.iter())
            .map(|(a, b)| a | b)
            .collect();
        FingerprintBits {
            words,
            nbits: self.nbits,
        }
    }

    /// Raw word slice for SIMD operations.
    pub fn words(&self) -> &[u64] {
        &self.words
    }

    pub fn nbits(&self) -> usize {
        self.nbits
    }

    /// Encode as lowercase hex string (for FPS format).
    pub fn to_hex(&self) -> String {
        let nbytes = self.nbits.div_ceil(8);
        let mut out = String::with_capacity(nbytes * 2);
        for i in 0..nbytes {
            let word = self.words[i / 8];
            let byte = (word >> ((i % 8) * 8)) as u8;
            out.push_str(&format!("{:02x}", byte));
        }
        out
    }

    /// Parse a hex string back into a FingerprintBits.
    pub fn from_hex(hex: &str, nbits: usize) -> Result<Self, String> {
        let mut fp = FingerprintBits::new(nbits);
        let nbytes = nbits.div_ceil(8);
        if hex.len() != nbytes * 2 {
            return Err(format!(
                "hex length {} doesn't match {} bits",
                hex.len(),
                nbits
            ));
        }
        for i in 0..nbytes {
            let byte = u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).map_err(|e| e.to_string())?;
            let word_idx = i / 8;
            let bit_shift = (i % 8) * 8;
            fp.words[word_idx] |= (byte as u64) << bit_shift;
        }
        Ok(fp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_get() {
        let mut fp = FingerprintBits::new(64);
        fp.set(0);
        fp.set(63);
        assert!(fp.get(0));
        assert!(fp.get(63));
        assert!(!fp.get(1));
    }

    #[test]
    fn test_count_ones() {
        let mut fp = FingerprintBits::new(128);
        fp.set(0);
        fp.set(64);
        fp.set(127);
        assert_eq!(fp.count_ones(), 3);
    }

    #[test]
    fn test_and_or() {
        let mut a = FingerprintBits::new(64);
        let mut b = FingerprintBits::new(64);
        a.set(0);
        a.set(1);
        b.set(1);
        b.set(2);

        let and = a.and(&b);
        assert!(and.get(1));
        assert!(!and.get(0));
        assert!(!and.get(2));

        let or = a.or(&b);
        assert!(or.get(0));
        assert!(or.get(1));
        assert!(or.get(2));
    }

    #[test]
    fn test_hex_roundtrip() {
        let mut fp = FingerprintBits::new(64);
        fp.set(0);
        fp.set(7);
        fp.set(8);
        let hex = fp.to_hex();
        let fp2 = FingerprintBits::from_hex(&hex, 64).unwrap();
        assert_eq!(fp, fp2);
    }

    #[test]
    fn test_folding_collision() {
        // Two different hashes mapping to same bit should still set it once
        let mut fp = FingerprintBits::new(16);
        fp.set(5);
        fp.set(5 + 16); // same bit after folding
        assert_eq!(fp.count_ones(), 1);
    }
}
