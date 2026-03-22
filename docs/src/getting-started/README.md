# Getting Started

molprint ships as both a CLI tool and a Rust library. Choose the approach that fits your workflow.

## Install the CLI

```bash
cargo install --git https://github.com/mariusrueve/molprint molprint-cli
```

This installs the `molprint` binary to `~/.cargo/bin/`. Make sure that directory is on your `PATH`.

## Use as a library

Add the relevant crates to your `Cargo.toml`:

```toml
[dependencies]
molprint-core = { git = "https://github.com/mariusrueve/molprint" }
molprint-fp   = { git = "https://github.com/mariusrueve/molprint" }
molprint-search = { git = "https://github.com/mariusrueve/molprint" }
```

You only need the crates relevant to your use case. If you just want to parse SMILES, `molprint-core` is sufficient. For fingerprints, add `molprint-fp`. For similarity search, add `molprint-search`.

## Build from source

```bash
git clone https://github.com/mariusrueve/molprint
cd molprint
cargo build --release
```

The `--release` flag is important for performance — debug builds can be 10–20× slower for fingerprint batch operations.

## Run tests

```bash
cargo test --workspace
```

## Run benchmarks

```bash
cargo bench
```

Benchmarks use [Criterion](https://bheisler.github.io/criterion.rs/book/) and output HTML reports to `target/criterion/`.
