use super::metrics::tanimoto;
use molprint_fp::FingerprintBits;
use rayon::prelude::*;

/// A single similarity search result.
#[derive(Debug, Clone)]
pub struct SearchHit {
    pub index: usize,
    pub similarity: f64,
}

/// Find all database entries with similarity ≥ threshold, sorted descending.
pub fn threshold_search(
    query: &FingerprintBits,
    database: &[FingerprintBits],
    threshold: f64,
) -> Vec<SearchHit> {
    let mut hits: Vec<SearchHit> = database
        .par_iter()
        .enumerate()
        .filter_map(|(i, fp)| {
            let sim = tanimoto(query, fp);
            if sim >= threshold {
                Some(SearchHit {
                    index: i,
                    similarity: sim,
                })
            } else {
                None
            }
        })
        .collect();

    hits.sort_unstable_by(|a, b| b.similarity.partial_cmp(&a.similarity).unwrap());
    hits
}

/// Find top-K most similar fingerprints.
pub fn top_k_search(
    query: &FingerprintBits,
    database: &[FingerprintBits],
    k: usize,
) -> Vec<SearchHit> {
    let mut hits: Vec<SearchHit> = database
        .par_iter()
        .enumerate()
        .map(|(i, fp)| SearchHit {
            index: i,
            similarity: tanimoto(query, fp),
        })
        .collect();

    hits.sort_unstable_by(|a, b| b.similarity.partial_cmp(&a.similarity).unwrap());
    hits.truncate(k);
    hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use molprint_core::smiles::parse_smiles;
    use molprint_fp::morgan::Morgan;
    use molprint_fp::traits::Fingerprinter;
    use molprint_fp::FingerprintBits;

    #[test]
    fn test_search_self_tanimoto() {
        let mol = parse_smiles("CCO").unwrap();
        let morgan = Morgan::new(2, 2048);
        let fp = morgan.fingerprint(&mol);
        let db = vec![fp.clone()];
        let hits = threshold_search(&fp, &db, 0.99);
        assert!(!hits.is_empty());
        assert!((hits[0].similarity - 1.0).abs() < 1e-9);
    }
}
