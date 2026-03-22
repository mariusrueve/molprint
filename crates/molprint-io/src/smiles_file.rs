use molprint_core::{smiles::parse_smiles, MolGraph};
use std::io::{BufRead, BufReader, Read};

/// Read a SMILES file (one per line, optional tab-separated ID).
/// Skips lines that fail to parse (logs warning to stderr).
pub struct SmilesFileReader<R: Read> {
    reader: BufReader<R>,
    line_number: usize,
}

impl<R: Read> SmilesFileReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            line_number: 0,
        }
    }
}

impl<R: Read> Iterator for SmilesFileReader<R> {
    type Item = (String, MolGraph);

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => return None,
                Err(_) => return None,
                Ok(_) => {}
            }
            self.line_number += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let (smiles, id) = if let Some(tab) = trimmed.find('\t') {
                (&trimmed[..tab], trimmed[tab + 1..].to_string())
            } else {
                (trimmed, format!("mol{}", self.line_number))
            };

            match parse_smiles(smiles) {
                Ok(mol) => return Some((id, mol)),
                Err(e) => {
                    eprintln!(
                        "Warning: line {}: failed to parse '{}': {}",
                        self.line_number, smiles, e
                    );
                    continue;
                }
            }
        }
    }
}
