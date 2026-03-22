use molprint_fp::FingerprintBits;
use std::io::{BufRead, BufReader, Read, Write};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FpsError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed FPS line: {0}")]
    BadLine(String),
    #[error("hex parse error: {0}")]
    HexError(String),
}

/// Write fingerprints to a chemfp-compatible FPS file.
pub struct FpsWriter<W: Write> {
    writer: W,
    #[allow(dead_code)]
    nbits: usize,
    #[allow(dead_code)]
    fp_type: String,
}

impl<W: Write> FpsWriter<W> {
    pub fn new(mut writer: W, nbits: usize, fp_type: &str) -> std::io::Result<Self> {
        writeln!(writer, "#FPS1")?;
        writeln!(writer, "#type={} nbits={}", fp_type, nbits)?;
        writeln!(writer, "#software=molprint/0.1.0")?;
        Ok(Self {
            writer,
            nbits,
            fp_type: fp_type.to_string(),
        })
    }

    /// Write an additional `#key=value` metadata line to the header.
    /// Must be called before any `write_fingerprint` calls.
    pub fn write_meta(&mut self, key: &str, value: &str) -> std::io::Result<()> {
        writeln!(self.writer, "#{}={}", key, value)
    }

    pub fn write_fingerprint(&mut self, id: &str, fp: &FingerprintBits) -> std::io::Result<()> {
        writeln!(self.writer, "{}\t{}", fp.to_hex(), id)
    }
}

/// Read fingerprints from a chemfp-compatible FPS file.
pub struct FpsReader<R: Read> {
    reader: BufReader<R>,
    nbits: usize,
    fp_type: String,
    meta: std::collections::HashMap<String, String>,
    /// First data line consumed during header parsing (if any).
    pending: Option<String>,
    header_read: bool,
}

impl<R: Read> FpsReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            nbits: 0,
            fp_type: String::new(),
            meta: std::collections::HashMap::new(),
            pending: None,
            header_read: false,
        }
    }

    /// The fingerprint type string from the `#type=` header line.
    pub fn fp_type(&self) -> &str {
        &self.fp_type
    }

    /// A metadata value from the header (e.g. `"radius"` → `"2"`).
    pub fn meta(&self, key: &str) -> Option<&str> {
        self.meta.get(key).map(String::as_str)
    }

    /// Read and process header, returning nbits.
    pub fn read_header(&mut self) -> Result<usize, FpsError> {
        if self.header_read {
            return Ok(self.nbits);
        }
        self.header_read = true;
        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(0) => return Ok(self.nbits),
                Err(e) => return Err(FpsError::Io(e)),
                Ok(_) => {}
            }
            let trimmed = line.trim();
            if !trimmed.starts_with('#') {
                // Save the first data line so the iterator can use it
                if !trimmed.is_empty() {
                    self.pending = Some(line);
                }
                break;
            }
            // Parse `#key=value` header fields; a line may have multiple
            // space-separated pairs, e.g. `#type=Morgan nbits=2048`
            let content = &trimmed[1..]; // strip leading '#'
            for token in content.split_whitespace() {
                if let Some(eq) = token.find('=') {
                    let key = token[..eq].to_string();
                    let value = token[eq + 1..].to_string();
                    match key.as_str() {
                        "type" => self.fp_type = value.clone(),
                        "nbits" | "num_bits" => {
                            if let Ok(n) = value.parse::<usize>() {
                                self.nbits = n;
                            }
                        }
                        _ => {}
                    }
                    self.meta.insert(key, value);
                }
            }
        }
        Ok(self.nbits)
    }

    fn parse_line(&self, line: &str) -> Option<Result<(String, FingerprintBits), FpsError>> {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            return None;
        }
        let parts: Vec<&str> = trimmed.splitn(2, '\t').collect();
        if parts.len() < 2 {
            return Some(Err(FpsError::BadLine(trimmed.to_string())));
        }
        let hex = parts[0];
        let id = parts[1].to_string();
        let nbits = if self.nbits > 0 {
            self.nbits
        } else {
            hex.len() * 4
        };
        Some(match FingerprintBits::from_hex(hex, nbits) {
            Ok(fp) => Ok((id, fp)),
            Err(e) => Err(FpsError::HexError(e)),
        })
    }
}

impl<R: Read> Iterator for FpsReader<R> {
    type Item = Result<(String, FingerprintBits), FpsError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Auto-read header on first call if not done yet
        if !self.header_read {
            if let Err(e) = self.read_header() {
                return Some(Err(e));
            }
        }

        // Serve pending line saved during header read
        if let Some(line) = self.pending.take() {
            if let Some(result) = self.parse_line(&line) {
                return Some(result);
            }
        }

        let mut line = String::new();
        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => return None,
                Err(e) => return Some(Err(FpsError::Io(e))),
                Ok(_) => {}
            }
            if let Some(result) = self.parse_line(&line) {
                return Some(result);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molprint_fp::FingerprintBits;

    #[test]
    fn test_fps_roundtrip() {
        let mut fp = FingerprintBits::new(64);
        fp.set(0);
        fp.set(7);
        fp.set(63);
        let hex = fp.to_hex();
        let fp2 = FingerprintBits::from_hex(&hex, 64).unwrap();
        assert_eq!(fp, fp2);

        // Test through FpsWriter/FpsReader
        let mut buf = Vec::new();
        let mut writer = FpsWriter::new(&mut buf, 64, "test").unwrap();
        writer.write_fingerprint("mol1", &fp).unwrap();
        drop(writer);

        let mut reader = FpsReader::new(buf.as_slice());
        let nbits = reader.read_header().unwrap();
        assert_eq!(nbits, 64);
        let (id, fp_read) = reader.next().unwrap().unwrap();
        assert_eq!(id, "mol1");
        assert_eq!(fp, fp_read);
    }
}

#[cfg(test)]
mod integration_tests {
    use super::*;
    use molprint_core::smiles::parse_smiles;
    use molprint_fp::morgan::Morgan;
    use molprint_fp::traits::Fingerprinter;

    #[test]
    fn test_fps_roundtrip_2048() {
        let mol = parse_smiles("CCO").unwrap();
        let morgan = Morgan::new(2, 2048);
        let fp_orig = morgan.fingerprint(&mol);

        let mut buf = Vec::new();
        let mut writer = FpsWriter::new(&mut buf, 2048, "Morgan").unwrap();
        writer.write_fingerprint("ethanol", &fp_orig).unwrap();
        drop(writer);

        let mut reader = FpsReader::new(buf.as_slice());
        let nbits = reader.read_header().unwrap();
        assert_eq!(nbits, 2048);
        let (id, fp_read) = reader.next().unwrap().unwrap();
        assert_eq!(id, "ethanol");
        assert_eq!(fp_orig.count_ones(), fp_read.count_ones());
        assert_eq!(fp_orig, fp_read, "fingerprint changed after FPS round-trip");
    }
}
