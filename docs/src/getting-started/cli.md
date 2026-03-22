# CLI Usage

The `molprint` CLI has two subcommands: `fp` for computing fingerprints, and `search` for similarity search.

## `molprint fp` — compute fingerprints

Reads a molecule file, computes fingerprints, and writes a chemfp-compatible FPS file.

```bash
molprint fp --input molecules.smi --output molecules.fps
```

### Options

| Flag | Default | Description |
|---|---|---|
| `--input` / `-i` | required | Input file (`.smi`, `.sdf`, `.sdf.gz`) |
| `--output` / `-o` | required | Output FPS file path |
| `--fp-type` / `-t` | `morgan` | Fingerprint type: `morgan` or `maccs` |
| `--radius` / `-r` | `2` | Morgan radius (ignored for MACCS) |
| `--nbits` / `-n` | `2048` | Bit width for Morgan (ignored for MACCS) |

### Examples

```bash
# Default: Morgan ECFP4, 2048 bits
molprint fp --input molecules.smi --output molecules.fps

# MACCS-166 structural keys
molprint fp --input molecules.smi --fp-type maccs --output maccs.fps

# Morgan ECFP6, 4096 bits, from SDF
molprint fp --input molecules.sdf --radius 3 --nbits 4096 --output ecfp6.fps

# Gzip-compressed SDF (auto-detected by extension)
molprint fp --input chembl.sdf.gz --output chembl.fps
```

### Input formats

File format is detected from the extension:

- `.smi` — SMILES file, one molecule per line. Format: `SMILES<tab>ID` or bare `SMILES`. If no ID is given, a sequential `mol1`, `mol2`, ... ID is assigned.
- `.sdf` — MDL V2000 SDF. The molecule name (first line of the mol block) is used as the ID.
- `.sdf.gz` — gzip-compressed SDF, read transparently via flate2.

### Output format

The output is a [chemfp](https://chemfp.com/)-compatible FPS v1 file:

```
#FPS1
#type=Morgan nbits=2048
#radius=2
#software=molprint/0.1.0
<hex_fingerprint>\t<id>
<hex_fingerprint>\t<id>
...
```

Fingerprints are encoded as lowercase hexadecimal, least-significant byte first (chemfp convention).

### Progress output

The CLI prints a progress summary to stderr when done:

```
Processed 10000 molecules in 0.14s (71428 mol/s)
```

## `molprint search` — similarity search

Searches a fingerprint database for molecules similar to a query.

```bash
molprint search --query "c1ccccc1" --db molecules.fps --threshold 0.4
```

### Options

| Flag | Default | Description |
|---|---|---|
| `--query` / `-q` | required | Query as a SMILES string |
| `--db` / `-d` | required | FPS database file path |
| `--threshold` | `0.7` | Tanimoto threshold for filtering |
| `--top-k` / `-k` | `10` | Return top-k results |

### How threshold and top-k interact

- If `--top-k` is non-zero (the default is 10), a top-k search is performed and then results below `--threshold` are filtered out.
- To get all results above a threshold without a hard top-k cap, set `--top-k 0`.

### Examples

```bash
# All compounds with Tanimoto >= 0.4 to benzene (no top-k cap)
molprint search --query "c1ccccc1" --db molecules.fps --threshold 0.4 --top-k 0

# Top 20 most similar, filtered to Tanimoto >= 0.7
molprint search --query "c1ccccc1" --db molecules.fps --top-k 20

# Top 5 most similar (no threshold filter)
molprint search --query "c1ccccc1" --db molecules.fps --top-k 5 --threshold 0.0
```

### Output format

Tab-separated `id<tab>similarity`, sorted by descending score:

```
CHEMBL123    0.8500
CHEMBL456    0.7800
CHEMBL789    0.7200
```

The fingerprint type (morgan/maccs) and parameters (radius, nbits) are read from the FPS header automatically, so the query fingerprint is computed with the same settings as the database.
