# Library Usage

This page walks through using molprint as a Rust library, from parsing a single SMILES string to running a parallel similarity search over thousands of compounds.

## Parsing SMILES

```rust
use molprint_core::smiles::parse_smiles;

let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
println!("{} atoms, {} bonds", mol.node_count(), mol.edge_count());
// 13 atoms, 13 bonds
```

`parse_smiles` returns a `Result<MolGraph, ParseError>`. In library code, propagate the error rather than unwrapping.

## Inspecting atoms

```rust
use molprint_core::smiles::parse_smiles;
use molprint_core::mol::graph::MolGraphExt;
use petgraph::graph::NodeIndex;

let mol = parse_smiles("CCO").unwrap(); // ethanol

for i in 0..mol.node_count() {
    let idx = NodeIndex::new(i);
    let atom = mol.atom(idx);
    println!(
        "atom {}: {:?} charge={} h_count={} in_ring={}",
        i, atom.element, atom.charge, atom.h_count, atom.in_ring
    );
}
```

## Computing fingerprints

### Morgan (ECFP)

```rust
use molprint_core::smiles::parse_smiles;
use molprint_fp::{morgan::Morgan, traits::Fingerprinter};

let mol = parse_smiles("c1ccccc1").unwrap();

// ECFP4 (radius=2, 2048 bits)
let morgan = Morgan::new(2, 2048);
let fp = morgan.fingerprint(&mol);

println!("{} bits set out of 2048", fp.count_ones());
println!("hex: {}", fp.to_hex());
```

`Morgan::new(radius, nbits)` where:
- `radius = 1` → ECFP2
- `radius = 2` → ECFP4 (default)
- `radius = 3` → ECFP6

### MACCS-166

```rust
use molprint_core::smiles::parse_smiles;
use molprint_fp::{maccs::Maccs166, traits::Fingerprinter};

let mol = parse_smiles("c1ccccc1").unwrap();
let maccs = Maccs166::new();
let fp = maccs.fingerprint(&mol);

println!("{} bits set out of 166", fp.count_ones());
```

MACCS-166 always produces 167-bit fingerprints (bit 0 is always 0; bits 1–166 correspond to the 166 structural keys).

### The `Fingerprinter` trait

Both `Morgan` and `Maccs166` implement `Fingerprinter`:

```rust
pub trait Fingerprinter {
    fn fingerprint(&self, mol: &MolGraph) -> FingerprintBits;
    fn name(&self) -> &str;
    fn nbits(&self) -> usize;
}
```

You can use `Box<dyn Fingerprinter>` to select the algorithm at runtime.

## Computing similarity

```rust
use molprint_core::smiles::parse_smiles;
use molprint_fp::{morgan::Morgan, traits::Fingerprinter};
use molprint_search::metrics::{tanimoto, dice, cosine};

let benzene = parse_smiles("c1ccccc1").unwrap();
let pyridine = parse_smiles("c1ccncc1").unwrap();

let morgan = Morgan::new(2, 2048);
let fp_a = morgan.fingerprint(&benzene);
let fp_b = morgan.fingerprint(&pyridine);

println!("Tanimoto: {:.3}", tanimoto(&fp_a, &fp_b));
println!("Dice:     {:.3}", dice(&fp_a, &fp_b));
println!("Cosine:   {:.3}", cosine(&fp_a, &fp_b));
```

All three metrics return `f64` in `[0.0, 1.0]`. Tanimoto is the standard in cheminformatics.

## Parallel screening

```rust
use molprint_core::smiles::parse_smiles;
use molprint_fp::{morgan::Morgan, traits::Fingerprinter};
use molprint_search::screen::{threshold_search, top_k_search};

let morgan = Morgan::new(2, 2048);

// Build a database (in practice, loaded from an FPS file)
let smiles_list = vec!["CCO", "CCCO", "c1ccccc1", "CC(=O)O", "c1ccncc1"];
let db_fps: Vec<_> = smiles_list.iter()
    .map(|smi| morgan.fingerprint(&parse_smiles(smi).unwrap()))
    .collect();

// Query
let query = parse_smiles("c1ccccc1").unwrap();
let query_fp = morgan.fingerprint(&query);

// Threshold search: all compounds with Tanimoto >= 0.5
let hits = threshold_search(&query_fp, &db_fps, 0.5);
for hit in &hits {
    println!("index={} similarity={:.3}", hit.index, hit.similarity);
}

// Top-k search: 3 most similar
let top3 = top_k_search(&query_fp, &db_fps, 3);
```

Both functions use Rayon for parallelism. They return `Vec<SearchHit>` sorted by descending similarity.

## Reading and writing FPS files

```rust
use std::fs::File;
use std::io::BufWriter;
use molprint_io::fps::{FpsWriter, FpsReader};

// Write
let out = File::create("output.fps").unwrap();
let mut writer = FpsWriter::new(BufWriter::new(out), 2048, "Morgan").unwrap();
writer.write_fingerprint("mol1", &fp).unwrap();

// Read
let f = File::open("output.fps").unwrap();
let mut reader = FpsReader::new(f);
let nbits = reader.read_header().unwrap();

for item in &mut reader {
    let (id, fp) = item.unwrap();
    println!("{}: {} bits set", id, fp.count_ones());
}
```

## Reading SMILES files

```rust
use std::fs::File;
use molprint_io::smiles_file::SmilesFileReader;

let f = File::open("molecules.smi").unwrap();
for (id, mol) in SmilesFileReader::new(f) {
    println!("{}: {} atoms", id, mol.node_count());
}
```

## Reading SDF files

```rust
use molprint_io::sdf::open_sdf;

for result in open_sdf("molecules.sdf").unwrap() {
    let record = result.unwrap();
    println!("{}: {} atoms", record.name, record.mol.node_count());
    // record.properties: HashMap<String, String>
}
```

`open_sdf` automatically handles `.sdf.gz` based on the file extension.

## SMARTS substructure matching

```rust
use molprint_core::smarts;
use molprint_core::smiles::parse_smiles;

let carbonyl = smarts::compile("C=O").unwrap();
let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();

if smarts::has_match(&aspirin, &carbonyl) {
    println!("aspirin contains a carbonyl group");
}

// Count all matches (returns one mapping per unique anchor atom)
let matches = smarts::subgraph_match(&aspirin, &carbonyl);
println!("{} carbonyl groups found", matches.len());
```
