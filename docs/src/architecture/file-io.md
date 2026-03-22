# File I/O

The `molprint-io` crate handles reading molecules from SMILES and SDF files, and reading/writing fingerprints in the FPS format.

## SMILES files (`smiles_file.rs`)

`SmilesFileReader` is an iterator over `(id, MolGraph)` pairs from a SMILES file:

```rust
use std::fs::File;
use molprint_io::smiles_file::SmilesFileReader;

let f = File::open("molecules.smi").unwrap();
for (id, mol) in SmilesFileReader::new(f) {
    println!("{}: {} atoms", id, mol.node_count());
}
```

### Input format

Accepts two formats per line:

- `SMILES<tab>ID` — SMILES and ID separated by a tab
- `SMILES` — bare SMILES; a sequential ID (`mol1`, `mol2`, …) is assigned

Lines starting with `#` are skipped as comments. Empty lines are skipped.

### Error handling

Lines where the SMILES fails to parse are silently skipped in the iterator (consistent with how large SMILES files from external sources are handled — one bad record shouldn't abort processing 10 million molecules). The CLI logs skipped records to stderr.

## SDF files (`sdf.rs`)

`open_sdf` returns a streaming iterator over SDF records:

```rust
use molprint_io::sdf::open_sdf;

for result in open_sdf("molecules.sdf").unwrap() {
    match result {
        Ok(record) => {
            println!("{}: {} atoms", record.name, record.mol.node_count());
            if let Some(mw) = record.properties.get("MW") {
                println!("  MW = {}", mw);
            }
        }
        Err(e) => eprintln!("Skipping record: {}", e),
    }
}
```

### SDF record

```rust
pub struct SdfRecord {
    pub name: String,                          // first line of the mol block
    pub mol: MolGraph,                         // parsed molecular graph
    pub properties: HashMap<String, String>,   // > <key> / value pairs
}
```

### Gzip support

`open_sdf` detects `.sdf.gz` by file extension and transparently wraps the file in a `flate2::read::GzDecoder`. No separate API is needed:

```rust
// Both work the same way
open_sdf("molecules.sdf")
open_sdf("molecules.sdf.gz")
```

### V2000 parsing

The parser handles MDL V2000 format mol blocks. The atom table and bond table are read line by line. Atom coordinates (x, y, z) are parsed but discarded — molprint is a 2D/topological library and doesn't use 3D coordinates.

Stereo bond wedge/dash indicators in the bond table are parsed and ignored. Properties outside the mol block (the `> <name>` / value pairs before `$$$$`) are collected into the `properties` map.

## FPS files (`fps.rs`)

The FPS format is chemfp's standard fingerprint exchange format. molprint implements a subset that is read/write compatible with chemfp.

### Writing

```rust
use molprint_io::fps::FpsWriter;
use std::io::BufWriter;

let out = File::create("output.fps").unwrap();
let mut writer = FpsWriter::new(BufWriter::new(out), 2048, "Morgan").unwrap();
writer.write_meta("radius", "2").unwrap();
writer.write_fingerprint("mol1", &fp).unwrap();
```

The writer outputs the FPS v1 header on construction:

```
#FPS1
#type=Morgan nbits=2048
#software=molprint/0.1.0
```

Additional metadata lines can be written with `write_meta` before any fingerprint records. The CLI uses this to store the Morgan radius so that `search` can reconstruct the exact fingerprinter from the database file.

Each fingerprint record is: `<hex>\t<id>\n`.

### Reading

```rust
use molprint_io::fps::FpsReader;

let mut reader = FpsReader::new(File::open("output.fps").unwrap());
let nbits = reader.read_header().unwrap(); // reads header, returns nbits

// Access metadata from header
let radius = reader.meta("radius"); // Option<&str>
let fp_type = reader.fp_type();     // &str

for item in &mut reader {
    let (id, fp) = item.unwrap();
}
```

`read_header` parses all `#key=value` comment lines at the top of the file. The reader handles the case where `nbits` is not in the header (infers from hex string length). If `read_header` is not called before iterating, it is called automatically on the first `next()` call.

### FPS format specification

```
#FPS1                                  ← magic header
#type=<name> nbits=<n>                 ← fingerprint type
#<key>=<value>                         ← optional metadata
<hex_fingerprint>\t<molecule_id>       ← one record per line
```

Hexadecimal encoding: bytes are written least-significant-first (little-endian), matching chemfp. For a 64-bit fingerprint with only bit 0 set, the hex is `0100000000000000` (8 bytes = 16 hex chars).
