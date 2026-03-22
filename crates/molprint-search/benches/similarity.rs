use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use molprint_fp::FingerprintBits;
use molprint_search::metrics::{cosine, dice, tanimoto};
use molprint_search::screen::{threshold_search, top_k_search};

fn make_fp(nbits: usize, step: usize) -> FingerprintBits {
    let mut fp = FingerprintBits::new(nbits);
    for i in (0..nbits).step_by(step) {
        fp.set(i);
    }
    fp
}

fn bench_metrics(c: &mut Criterion) {
    let mut group = c.benchmark_group("similarity_metrics");

    for nbits in [512usize, 1024, 2048, 4096] {
        let a = make_fp(nbits, 3);
        let b = make_fp(nbits, 5);

        group.throughput(Throughput::Elements(nbits as u64));

        group.bench_with_input(BenchmarkId::new("tanimoto", nbits), &nbits, |bench, _| {
            bench.iter(|| tanimoto(&a, &b))
        });
        group.bench_with_input(BenchmarkId::new("dice", nbits), &nbits, |bench, _| {
            bench.iter(|| dice(&a, &b))
        });
        group.bench_with_input(BenchmarkId::new("cosine", nbits), &nbits, |bench, _| {
            bench.iter(|| cosine(&a, &b))
        });
    }
    group.finish();
}

fn bench_screening(c: &mut Criterion) {
    let nbits = 2048;

    // Build a database of 10 000 fingerprints with varying density
    let db: Vec<FingerprintBits> = (1..=10_000)
        .map(|i| make_fp(nbits, (i % 7) + 2)) // step 2-8, varied density
        .collect();

    let query = make_fp(nbits, 3);

    let mut group = c.benchmark_group("screening");
    group.throughput(Throughput::Elements(db.len() as u64));

    group.bench_function("threshold_0.5_10k", |b| {
        b.iter(|| threshold_search(&query, &db, 0.5))
    });

    group.bench_function("threshold_0.3_10k", |b| {
        b.iter(|| threshold_search(&query, &db, 0.3))
    });

    group.bench_function("top_k_10_10k", |b| b.iter(|| top_k_search(&query, &db, 10)));

    group.bench_function("top_k_100_10k", |b| {
        b.iter(|| top_k_search(&query, &db, 100))
    });

    group.finish();
}

fn bench_screening_scale(c: &mut Criterion) {
    let nbits = 2048;
    let query = make_fp(nbits, 3);
    let mut group = c.benchmark_group("screening_scale");

    for size in [1_000usize, 10_000, 100_000] {
        let db: Vec<FingerprintBits> = (1..=size).map(|i| make_fp(nbits, (i % 7) + 2)).collect();

        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(BenchmarkId::new("threshold_0.4", size), &size, |b, _| {
            b.iter(|| threshold_search(&query, &db, 0.4))
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_metrics,
    bench_screening,
    bench_screening_scale
);
criterion_main!(benches);
