use clap::Parser;
use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;

use molprint_core::smiles::parse_smiles;
use molprint_fp::maccs::Maccs166;
use molprint_fp::morgan::Morgan;
use molprint_fp::traits::Fingerprinter;
use molprint_io::fps::{FpsReader, FpsWriter};
use molprint_io::sdf::open_sdf;
use molprint_io::smiles_file::SmilesFileReader;
use molprint_search::screen::{threshold_search, top_k_search};

#[derive(Parser)]
#[command(
    name = "molprint",
    about = "Molecular fingerprints & similarity search"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(clap::Subcommand)]
enum Commands {
    /// Compute fingerprints for molecules
    Fp {
        #[arg(short, long)]
        input: String,
        #[arg(short = 't', long, default_value = "morgan")]
        fp_type: String,
        #[arg(short, long, default_value_t = 2)]
        radius: u32,
        #[arg(short, long, default_value_t = 2048)]
        nbits: usize,
        #[arg(short, long)]
        output: String,
    },
    /// Search a fingerprint database for similar molecules
    Search {
        #[arg(short, long)]
        query: String,
        #[arg(short, long)]
        db: String,
        #[arg(long, default_value_t = 0.7)]
        threshold: f64,
        #[arg(short, long, default_value_t = 10)]
        top_k: usize,
    },
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Fp {
            input,
            fp_type,
            radius,
            nbits,
            output,
        } => {
            let fingerprinter: Box<dyn Fingerprinter> = match fp_type.as_str() {
                "morgan" | "ecfp" => Box::new(Morgan::new(radius, nbits)),
                "maccs" => Box::new(Maccs166::new()),
                other => {
                    eprintln!(
                        "Unknown fingerprint type: {}. Use 'morgan' or 'maccs'.",
                        other
                    );
                    std::process::exit(1);
                }
            };

            let out_file = File::create(&output).unwrap_or_else(|e| {
                eprintln!("Cannot create output file {}: {}", output, e);
                std::process::exit(1);
            });

            let mut writer = FpsWriter::new(
                BufWriter::new(out_file),
                fingerprinter.nbits(),
                fingerprinter.name(),
            )
            .unwrap_or_else(|e| {
                eprintln!("Cannot write FPS header: {}", e);
                std::process::exit(1);
            });

            // Write radius metadata for Morgan so search can reconstruct the fingerprinter
            if fp_type == "morgan" || fp_type == "ecfp" {
                writer
                    .write_meta("radius", &radius.to_string())
                    .unwrap_or_else(|e| {
                        eprintln!("Cannot write FPS metadata: {}", e);
                        std::process::exit(1);
                    });
            }

            let reader_file = File::open(&input).unwrap_or_else(|e| {
                eprintln!("Cannot open input file {}: {}", input, e);
                std::process::exit(1);
            });

            let start = Instant::now();
            let mut count = 0usize;

            // Auto-detect SDF (plain or gzip) vs SMILES by file extension
            let lc = input.to_lowercase();
            let is_sdf = lc.ends_with(".sdf") || lc.ends_with(".sdf.gz");
            if is_sdf {
                let sdf_reader = open_sdf(&input).unwrap_or_else(|e| {
                    eprintln!("Cannot open SDF file {}: {}", input, e);
                    std::process::exit(1);
                });
                for result in sdf_reader {
                    match result {
                        Ok(record) => {
                            let id = if record.name.is_empty() {
                                format!("mol{}", count + 1)
                            } else {
                                record.name
                            };
                            let fp = fingerprinter.fingerprint(&record.mol);
                            writer.write_fingerprint(&id, &fp).unwrap_or_else(|e| {
                                eprintln!("Write error: {}", e);
                            });
                            count += 1;
                        }
                        Err(e) => {
                            eprintln!("Warning: skipping SDF record: {}", e);
                        }
                    }
                }
            } else {
                for (id, mol) in SmilesFileReader::new(reader_file) {
                    let fp = fingerprinter.fingerprint(&mol);
                    writer.write_fingerprint(&id, &fp).unwrap_or_else(|e| {
                        eprintln!("Write error: {}", e);
                    });
                    count += 1;
                }
            }

            let elapsed = start.elapsed();
            let throughput = if elapsed.as_secs_f64() > 0.0 {
                count as f64 / elapsed.as_secs_f64()
            } else {
                f64::INFINITY
            };
            eprintln!(
                "Processed {} molecules in {:.2}s ({:.0} mol/s)",
                count,
                elapsed.as_secs_f64(),
                throughput
            );
        }

        Commands::Search {
            query,
            db,
            threshold,
            top_k,
        } => {
            // Parse query SMILES
            let query_mol = parse_smiles(&query).unwrap_or_else(|e| {
                eprintln!("Cannot parse query SMILES '{}': {}", query, e);
                std::process::exit(1);
            });

            // Load fingerprint database
            let db_file = File::open(&db).unwrap_or_else(|e| {
                eprintln!("Cannot open database file {}: {}", db, e);
                std::process::exit(1);
            });

            let mut fps_reader = FpsReader::new(db_file);
            let nbits = fps_reader.read_header().unwrap_or(2048);
            let nbits = if nbits == 0 { 2048 } else { nbits };

            let mut db_fps = Vec::new();
            let mut db_ids = Vec::new();

            for item in &mut fps_reader {
                match item {
                    Ok((id, fp)) => {
                        db_ids.push(id);
                        db_fps.push(fp);
                    }
                    Err(e) => {
                        eprintln!("Warning: skipping FPS record: {}", e);
                    }
                }
            }

            if db_fps.is_empty() {
                eprintln!("No fingerprints loaded from {}", db);
                std::process::exit(1);
            }

            // Reconstruct fingerprinter matching the database
            let db_fp_type = fps_reader.fp_type().to_lowercase();
            let db_radius: u32 = fps_reader
                .meta("radius")
                .and_then(|v| v.parse().ok())
                .unwrap_or(2);

            let fingerprinter: Box<dyn Fingerprinter> = if db_fp_type.starts_with("maccs") {
                Box::new(Maccs166::new())
            } else {
                // Default to Morgan for "morgan", "ecfp", or unknown types
                Box::new(Morgan::new(db_radius, nbits))
            };

            let query_fp = fingerprinter.fingerprint(&query_mol);

            // Run search
            let hits = if top_k > 0 {
                top_k_search(&query_fp, &db_fps, top_k)
                    .into_iter()
                    .filter(|h| h.similarity >= threshold)
                    .collect::<Vec<_>>()
            } else {
                threshold_search(&query_fp, &db_fps, threshold)
            };

            for hit in &hits {
                println!("{}\t{:.4}", db_ids[hit.index], hit.similarity);
            }

            if hits.is_empty() {
                eprintln!("No hits found above threshold {}", threshold);
            }
        }
    }
}
