use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use molprint_core::smiles::parse_smiles;
use molprint_fp::maccs::Maccs166;
use molprint_fp::morgan::Morgan;
use molprint_fp::traits::Fingerprinter;

const SMILES: &[&str] = &[
    "C",
    "CCO",
    "c1ccccc1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "CC(=O)Nc1ccc(O)cc1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
    "CCN(CC)CC(=O)Nc1c(C)cccc1C",
    "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
    "CC(C)c1c(C(=O)Nc2ccccc2F)c(c2ccc(F)cc2)c(n1CCc1ccc(O)cc1)CC(O)CC(=O)O",
    "CCCc1nn(c2ccc(cc12)S(=O)(=O)N1CCN(CC1)C)c1ccncc1C(=O)OCC",
];

fn make_mols() -> Vec<molprint_core::MolGraph> {
    SMILES.iter().map(|s| parse_smiles(s).unwrap()).collect()
}

fn bench_morgan_radii(c: &mut Criterion) {
    let mols = make_mols();
    let mut group = c.benchmark_group("morgan_radius");

    for radius in [1u32, 2, 3] {
        let fp = Morgan::new(radius, 2048);
        group.bench_with_input(BenchmarkId::new("r", radius), &radius, |b, _| {
            b.iter(|| {
                for mol in &mols {
                    let _ = fp.fingerprint(mol);
                }
            })
        });
    }
    group.finish();
}

fn bench_morgan_nbits(c: &mut Criterion) {
    let mols = make_mols();
    let mut group = c.benchmark_group("morgan_nbits");

    for nbits in [512usize, 1024, 2048, 4096] {
        let fp = Morgan::new(2, nbits);
        group.bench_with_input(BenchmarkId::new("bits", nbits), &nbits, |b, _| {
            b.iter(|| {
                for mol in &mols {
                    let _ = fp.fingerprint(mol);
                }
            })
        });
    }
    group.finish();
}

fn bench_maccs(c: &mut Criterion) {
    let mols = make_mols();
    let fp = Maccs166::new();

    c.bench_function("maccs166_batch_12mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = fp.fingerprint(mol);
            }
        })
    });
}

fn bench_morgan_throughput(c: &mut Criterion) {
    // Sustained throughput: aspirin × 1000
    let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let morgan = Morgan::new(2, 2048);

    c.benchmark_group("morgan_throughput")
        .throughput(Throughput::Elements(1000))
        .bench_function("aspirin_ecfp4_x1000", |b| {
            b.iter(|| {
                for _ in 0..1000 {
                    let _ = morgan.fingerprint(&mol);
                }
            })
        });
}

fn bench_fingerprint_per_mol(c: &mut Criterion) {
    let mols = make_mols();
    let morgan = Morgan::new(2, 2048);
    let maccs = Maccs166::new();
    let mut group = c.benchmark_group("fingerprint_per_mol");

    for (i, (mol, smi)) in mols.iter().zip(SMILES.iter()).enumerate() {
        if i % 3 == 0 {
            // sample every 3rd to keep output manageable
            group.bench_with_input(
                BenchmarkId::new("morgan", &smi[..smi.len().min(20)]),
                mol,
                |b, mol| b.iter(|| morgan.fingerprint(mol)),
            );
            group.bench_with_input(
                BenchmarkId::new("maccs", &smi[..smi.len().min(20)]),
                mol,
                |b, mol| b.iter(|| maccs.fingerprint(mol)),
            );
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_morgan_radii,
    bench_morgan_nbits,
    bench_maccs,
    bench_morgan_throughput,
    bench_fingerprint_per_mol,
);
criterion_main!(benches);
