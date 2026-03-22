use flate2::read::GzDecoder;
use molprint_core::mol::graph::MolGraphExt;
use molprint_core::ring::assign_ring_info;
use molprint_core::{Atom, BondType, Element, MolGraph};
use petgraph::graph::NodeIndex;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SdfError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed counts line at record {0}")]
    BadCountsLine(usize),
    #[error("malformed atom line at record {0}, line {1}")]
    BadAtomLine(usize, usize),
    #[error("malformed bond line at record {0}, line {1}")]
    BadBondLine(usize, usize),
}

/// A parsed SDF record.
pub struct SdfRecord {
    pub name: String,
    pub mol: MolGraph,
    pub properties: std::collections::HashMap<String, String>,
}

/// Streaming SDF V2000 reader. Yields one `SdfRecord` per molecule.
pub struct SdfReader<R: Read> {
    reader: BufReader<R>,
    record_num: usize,
}

impl<R: Read> SdfReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            record_num: 0,
        }
    }

    fn read_record(&mut self) -> Option<Result<SdfRecord, SdfError>> {
        self.record_num += 1;
        let rec = self.record_num;

        // Line 1: molecule name
        let mut name_line = String::new();
        match self.reader.read_line(&mut name_line) {
            Ok(0) => return None,
            Err(e) => return Some(Err(SdfError::Io(e))),
            Ok(_) => {}
        }
        let name = name_line.trim().to_string();

        // Lines 2-3: program/comment (skip)
        for _ in 0..2 {
            let mut buf = String::new();
            if let Err(e) = self.reader.read_line(&mut buf) {
                return Some(Err(SdfError::Io(e)));
            }
        }

        // Line 4: counts line "aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv"
        let mut counts_line = String::new();
        if let Err(e) = self.reader.read_line(&mut counts_line) {
            return Some(Err(SdfError::Io(e)));
        }
        let counts_line = counts_line.trim();
        if counts_line.len() < 6 {
            return Some(Err(SdfError::BadCountsLine(rec)));
        }

        let n_atoms: usize = match counts_line[0..3].trim().parse() {
            Ok(n) => n,
            Err(_) => return Some(Err(SdfError::BadCountsLine(rec))),
        };
        let n_bonds: usize = match counts_line[3..6].trim().parse() {
            Ok(n) => n,
            Err(_) => return Some(Err(SdfError::BadCountsLine(rec))),
        };

        // Atom block
        let mut graph = MolGraph::new_undirected();
        for atom_line_num in 0..n_atoms {
            let mut line = String::new();
            if let Err(e) = self.reader.read_line(&mut line) {
                return Some(Err(SdfError::Io(e)));
            }
            let line = line.trim_end();
            if line.len() < 34 {
                return Some(Err(SdfError::BadAtomLine(rec, atom_line_num)));
            }
            // Columns 31-34: element symbol (1-indexed: 32-34)
            let sym = line[31..34].trim();
            let element = Element::from_symbol(sym).unwrap_or(Element::Unknown);
            // Columns 36-39: charge code
            let charge_code: i8 = line
                .get(36..39)
                .and_then(|s| s.trim().parse::<u8>().ok())
                .map(sdf_charge_code_to_charge)
                .unwrap_or(0);
            let mut atom = Atom::new(element);
            atom.charge = charge_code;
            graph.add_node(atom);
        }

        // Bond block
        for bond_line_num in 0..n_bonds {
            let mut line = String::new();
            if let Err(e) = self.reader.read_line(&mut line) {
                return Some(Err(SdfError::Io(e)));
            }
            let line = line.trim_end();
            if line.len() < 9 {
                return Some(Err(SdfError::BadBondLine(rec, bond_line_num)));
            }
            let a1: usize = match line[0..3].trim().parse::<usize>() {
                Ok(n) if n > 0 => n - 1,
                _ => return Some(Err(SdfError::BadBondLine(rec, bond_line_num))),
            };
            let a2: usize = match line[3..6].trim().parse::<usize>() {
                Ok(n) if n > 0 => n - 1,
                _ => return Some(Err(SdfError::BadBondLine(rec, bond_line_num))),
            };
            let btype: u8 = match line[6..9].trim().parse() {
                Ok(n) => n,
                _ => return Some(Err(SdfError::BadBondLine(rec, bond_line_num))),
            };
            let bond = match btype {
                1 => BondType::Single,
                2 => BondType::Double,
                3 => BondType::Triple,
                4 => BondType::Aromatic,
                _ => BondType::Single,
            };
            if a1 < n_atoms && a2 < n_atoms {
                graph.add_edge(NodeIndex::new(a1), NodeIndex::new(a2), bond);
            }
        }

        graph.assign_implicit_hydrogens();
        assign_ring_info(&mut graph);

        // Read properties until "$$$$"
        let mut properties = std::collections::HashMap::new();
        let mut current_prop: Option<String> = None;
        let mut prop_value = String::new();

        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(0) => break,
                Err(e) => return Some(Err(SdfError::Io(e))),
                Ok(_) => {}
            }
            let trimmed = line.trim();
            if trimmed == "$$$$" {
                if let Some(key) = current_prop.take() {
                    properties.insert(key, prop_value.trim().to_string());
                    prop_value.clear();
                }
                break;
            }
            if trimmed.starts_with("> <") {
                if let Some(key) = current_prop.take() {
                    properties.insert(key, prop_value.trim().to_string());
                    prop_value.clear();
                }
                if let Some(end) = trimmed.find('>') {
                    current_prop = Some(trimmed[3..end].to_string());
                }
            } else if current_prop.is_some() {
                if !prop_value.is_empty() {
                    prop_value.push('\n');
                }
                prop_value.push_str(trimmed);
            }
        }

        Some(Ok(SdfRecord {
            name,
            mol: graph,
            properties,
        }))
    }
}

fn sdf_charge_code_to_charge(code: u8) -> i8 {
    match code {
        0 => 0,
        1 => 3,
        2 => 2,
        3 => 1,
        4 => 0, // doublet radical — treat as 0
        5 => -1,
        6 => -2,
        7 => -3,
        _ => 0,
    }
}

impl<R: Read> Iterator for SdfReader<R> {
    type Item = Result<SdfRecord, SdfError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_record()
    }
}

/// A type-erased SDF reader that works for both plain `.sdf` and gzip `.sdf.gz` files.
pub type AnySdfReader = SdfReader<Box<dyn Read>>;

/// Open an SDF file for streaming, transparently decompressing `.sdf.gz` files.
///
/// # Errors
/// Returns an error if the file cannot be opened.
pub fn open_sdf<P: AsRef<Path>>(path: P) -> std::io::Result<AnySdfReader> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let reader: Box<dyn Read> = if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(SdfReader::new(reader))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    use tempfile::NamedTempFile;

    const ETHANOL_SDF: &str = "\
ethanol\n\
  molprint\n\
\n\
  3  2  0  0  0  0  0  0  0  0  0 V2000\n\
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0\n\
  2  3  1  0\n\
$$$$\n";

    #[test]
    fn test_sdf_gz_roundtrip() {
        // Write a gzip-compressed SDF to a temp file
        let mut tmp = NamedTempFile::with_suffix(".sdf.gz").unwrap();
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(ETHANOL_SDF.as_bytes()).unwrap();
        let compressed = gz.finish().unwrap();
        tmp.write_all(&compressed).unwrap();
        tmp.flush().unwrap();

        // Read it back via open_sdf
        let mut reader = open_sdf(tmp.path()).unwrap();
        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.name, "ethanol");
        assert_eq!(record.mol.node_count(), 3);
        assert!(reader.next().is_none());
    }
}
