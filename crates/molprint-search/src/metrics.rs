use molprint_fp::FingerprintBits;

/// Tanimoto coefficient (Jaccard index for bit vectors).
/// T(A,B) = |A∩B| / |A∪B| = popcount(A&B) / popcount(A|B)
pub fn tanimoto(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    let intersection = a.and(b).count_ones() as f64;
    let union = a.or(b).count_ones() as f64;
    if union == 0.0 {
        0.0
    } else {
        intersection / union
    }
}

/// Dice coefficient: 2*|A∩B| / (|A| + |B|)
pub fn dice(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    let intersection = a.and(b).count_ones() as f64;
    let sum = (a.count_ones() + b.count_ones()) as f64;
    if sum == 0.0 {
        0.0
    } else {
        2.0 * intersection / sum
    }
}

/// Cosine similarity: |A∩B| / sqrt(|A| * |B|)
pub fn cosine(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    let intersection = a.and(b).count_ones() as f64;
    let prod = (a.count_ones() as f64) * (b.count_ones() as f64);
    if prod == 0.0 {
        0.0
    } else {
        intersection / prod.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fp_from_bits(bits: &[usize], nbits: usize) -> FingerprintBits {
        let mut fp = FingerprintBits::new(nbits);
        for &b in bits {
            fp.set(b);
        }
        fp
    }

    #[test]
    fn test_tanimoto_self_similarity() {
        let fp = fp_from_bits(&[0, 1, 5, 10], 64);
        assert!((tanimoto(&fp, &fp) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_tanimoto_all_zeros() {
        let a = FingerprintBits::new(64);
        let b = FingerprintBits::new(64);
        assert_eq!(tanimoto(&a, &b), 0.0);
    }

    #[test]
    fn test_tanimoto_known_value() {
        // A = {0,1}, B = {1,2}: intersection=1, union=3 → T=1/3
        let a = fp_from_bits(&[0, 1], 64);
        let b = fp_from_bits(&[1, 2], 64);
        assert!((tanimoto(&a, &b) - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_tanimoto_symmetry() {
        let a = fp_from_bits(&[0, 3, 7], 64);
        let b = fp_from_bits(&[3, 5, 7, 9], 64);
        assert!((tanimoto(&a, &b) - tanimoto(&b, &a)).abs() < 1e-9);
    }

    #[test]
    fn test_dice() {
        let a = fp_from_bits(&[0, 1], 64);
        let b = fp_from_bits(&[1, 2], 64);
        // 2*1 / (2+2) = 0.5
        assert!((dice(&a, &b) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_cosine() {
        let a = fp_from_bits(&[0, 1], 64);
        let b = fp_from_bits(&[1, 2], 64);
        // 1 / sqrt(2*2) = 0.5
        assert!((cosine(&a, &b) - 0.5).abs() < 1e-9);
    }
}

#[cfg(test)]
mod integration_tests {
    use super::*;
    use molprint_core::smiles::parse_smiles;
    use molprint_fp::morgan::Morgan;
    use molprint_fp::traits::Fingerprinter;

    #[test]
    fn test_ethanol_propanol_similarity() {
        let mol1 = parse_smiles("CCO").unwrap();
        let mol2 = parse_smiles("CCCO").unwrap();
        let morgan = Morgan::new(2, 2048);
        let fp1 = morgan.fingerprint(&mol1);
        let fp2 = morgan.fingerprint(&mol2);
        let sim = tanimoto(&fp1, &fp2);
        println!("CCO vs CCCO Tanimoto: {:.4}", sim);
        println!("CCO bits: {}", fp1.count_ones());
        println!("CCCO bits: {}", fp2.count_ones());
        assert!(sim > 0.0, "ethanol and propanol should share some bits");
    }
}
