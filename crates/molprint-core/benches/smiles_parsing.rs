use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use molprint_core::smiles::parse_smiles;

/// Representative molecules spanning a range of sizes and complexity.
const SMILES: &[(&str, &str)] = &[
    ("methane", "C"),
    ("ethanol", "CCO"),
    ("benzene", "c1ccccc1"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("caffeine", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
    ("lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    (
        "ciprofloxacin",
        "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
    ),
    (
        "atorvastatin",
        "CC(C)c1c(C(=O)Nc2ccccc2F)c(c2ccc(F)cc2)c(n1CCc1ccc(O)cc1)CC(O)CC(=O)O",
    ),
    (
        "sildenafil",
        "CCCc1nn(c2ccc(cc12)S(=O)(=O)N1CCN(CC1)C)c1ccncc1C(=O)OCC",
    ),
    (
        "taxol_core",
        "OC(C(=O)OC1C(OC(=O)c2ccccc2)C(C)(O)CC3(O)C(CC(=O)C(C(C)(C)O)C13)OC(=O)C)c1ccccc1",
    ),
    (
        "erythromycin_core",
        "CCC1OC(=O)C(CC(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(=O)C(C)C(O)C(C)(OC3OC(C)C(O)(C)C3OC)C1C)(C)O",
    ),
];

fn bench_parse_individual(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_smiles");

    for (name, smi) in SMILES {
        let heavy_atoms = parse_smiles(smi).map(|m| m.node_count()).unwrap_or(0);
        group.throughput(Throughput::Elements(heavy_atoms as u64));
        group.bench_with_input(BenchmarkId::new("atoms", name), smi, |b, smi| {
            b.iter(|| parse_smiles(smi).unwrap())
        });
    }
    group.finish();
}

fn bench_parse_batch(c: &mut Criterion) {
    let smiles: Vec<&str> = SMILES.iter().map(|(_, s)| *s).collect();

    c.bench_function("parse_smiles_batch_14mol", |b| {
        b.iter(|| {
            for smi in &smiles {
                let _ = parse_smiles(smi).unwrap();
            }
        })
    });
}

fn bench_parse_throughput(c: &mut Criterion) {
    // Repeat a mid-size molecule 1000× to measure sustained throughput
    let smi = "CC(=O)Oc1ccccc1C(=O)O"; // aspirin
    c.benchmark_group("parse_throughput")
        .throughput(Throughput::Elements(1000))
        .bench_function("aspirin_x1000", |b| {
            b.iter(|| {
                for _ in 0..1000 {
                    let _ = parse_smiles(smi).unwrap();
                }
            })
        });
}

criterion_group!(
    benches,
    bench_parse_individual,
    bench_parse_batch,
    bench_parse_throughput
);
criterion_main!(benches);
