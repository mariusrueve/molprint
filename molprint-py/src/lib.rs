use molprint_fp::bitvec::FingerprintBits;
use molprint_fp::maccs::Maccs166;
use molprint_fp::morgan::Morgan;
use molprint_fp::traits::Fingerprinter;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Convert a FingerprintBits to a Python bytes object (little-endian, byte-per-byte).
fn fp_to_bytes(fp: &FingerprintBits) -> Vec<u8> {
    let nbytes = fp.nbits().div_ceil(8);
    let words = fp.words();
    (0..nbytes)
        .map(|i| (words[i / 8] >> ((i % 8) * 8)) as u8)
        .collect()
}

/// Compute Tanimoto directly on raw byte slices (avoids hex round-trip).
fn tanimoto_bytes(a: &[u8], b: &[u8]) -> f64 {
    let mut and_count = 0u32;
    let mut or_count = 0u32;
    for (x, y) in a.iter().zip(b.iter()) {
        and_count += (x & y).count_ones();
        or_count += (x | y).count_ones();
    }
    if or_count == 0 {
        1.0
    } else {
        and_count as f64 / or_count as f64
    }
}

/// A parsed molecule.  Construct with `Mol.from_smiles(smiles)`.
#[pyclass(name = "Mol")]
struct PyMol {
    graph: molprint_core::MolGraph,
}

#[pymethods]
impl PyMol {
    /// Parse a SMILES string.  Raises `ValueError` on invalid input.
    #[staticmethod]
    fn from_smiles(smiles: &str) -> PyResult<Self> {
        molprint_core::smiles::parse_smiles(smiles)
            .map(|graph| PyMol { graph })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Compute the 166-bit MACCS fingerprint.
    ///
    /// Returns 21 bytes (166 bits, little-endian bit order).
    fn maccs(&self) -> Vec<u8> {
        fp_to_bytes(&Maccs166::new().fingerprint(&self.graph))
    }

    /// Compute a Morgan (ECFP-style) fingerprint.
    ///
    /// Args:
    ///     radius: neighbourhood radius (default 2)
    ///     nbits:  bit vector length (default 2048)
    ///
    /// Returns `nbits / 8` bytes (little-endian bit order).
    #[pyo3(signature = (radius=2, nbits=2048))]
    fn morgan(&self, radius: usize, nbits: usize) -> Vec<u8> {
        fp_to_bytes(&Morgan::new(radius, nbits).fingerprint(&self.graph))
    }
}

/// Tanimoto (Jaccard) similarity between two fingerprint byte strings.
///
/// Both byte strings must have the same length.
/// Returns a float in [0.0, 1.0].
#[pyfunction]
fn tanimoto(a: &[u8], b: &[u8]) -> PyResult<f64> {
    if a.len() != b.len() {
        return Err(PyValueError::new_err(format!(
            "fingerprint lengths differ: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    Ok(tanimoto_bytes(a, b))
}

/// Compute MACCS fingerprints for a list of SMILES strings (batch).
///
/// Invalid SMILES are silently skipped (None in output).
/// Returns a list of `bytes | None`, one per input SMILES.
#[pyfunction]
fn batch_maccs(smiles_list: Vec<String>) -> Vec<Option<Vec<u8>>> {
    let algo = Maccs166::new();
    smiles_list
        .iter()
        .map(|smi| {
            molprint_core::smiles::parse_smiles(smi)
                .ok()
                .map(|mol| fp_to_bytes(&algo.fingerprint(&mol)))
        })
        .collect()
}

/// Compute Morgan fingerprints for a list of SMILES strings (batch).
///
/// Invalid SMILES are silently skipped (None in output).
/// Returns a list of `bytes | None`, one per input SMILES.
#[pyfunction]
#[pyo3(signature = (smiles_list, radius=2, nbits=2048))]
fn batch_morgan(smiles_list: Vec<String>, radius: usize, nbits: usize) -> Vec<Option<Vec<u8>>> {
    let algo = Morgan::new(radius, nbits);
    smiles_list
        .iter()
        .map(|smi| {
            molprint_core::smiles::parse_smiles(smi)
                .ok()
                .map(|mol| fp_to_bytes(&algo.fingerprint(&mol)))
        })
        .collect()
}

/// molprint — high-performance molecular fingerprints.
///
/// Quick start::
///
///     import molprint
///
///     mol = molprint.Mol.from_smiles("c1ccccc1")
///     fp  = mol.maccs()                    # bytes, 21 bytes = 166 bits
///     fp2 = mol.morgan(radius=2, nbits=2048)
///
///     sim = molprint.tanimoto(fp, fp)      # 1.0
///
///     fps = molprint.batch_maccs(["CCO", "c1ccccc1"])   # list[bytes]
#[pymodule]
fn molprint(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMol>()?;
    m.add_function(wrap_pyfunction!(tanimoto, m)?)?;
    m.add_function(wrap_pyfunction!(batch_maccs, m)?)?;
    m.add_function(wrap_pyfunction!(batch_morgan, m)?)?;
    Ok(())
}
