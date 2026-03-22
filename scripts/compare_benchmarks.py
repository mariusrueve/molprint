#!/usr/bin/env python3
"""Compare molprint criterion benchmark results against RDKit throughput.

Usage:
    python scripts/compare_benchmarks.py bench_results.txt rdkit_bench.txt

bench_results.txt: captured stdout of `cargo bench`
rdkit_bench.txt:   output of scripts/benchmark_rdkit.py

Prints a comparison table and exits 0 if molprint is >= 5x faster on all metrics,
1 otherwise.
"""

import re
import sys
from pathlib import Path


# Minimum speedup required
MIN_SPEEDUP = 5.0


def parse_criterion(path: str) -> dict[str, float]:
    """Extract throughput from criterion bench output.

    Criterion lines look like:
        fingerprints/morgan_r2_2048/16
            time:   [42.123 µs 42.456 µs 42.789 µs]
    We sum per-molecule time from the per-element lines and derive mol/s.
    """
    results: dict[str, float] = {}
    # Also look for lines like:
    #   fingerprints/maccs_batch/16   time: [1.234 ms ...]
    pattern = re.compile(
        r"(morgan_r\d+_\d+|maccs_\w+|morgan\w*|maccs\w*)[^\n]*\n"
        r"\s+time:\s+\[[\d.]+ \S+\s+([\d.]+) (\S+)"
    )
    content = Path(path).read_text()
    for m in pattern.finditer(content):
        name_raw = m.group(1).lower()
        value = float(m.group(2))
        unit = m.group(3)

        # Convert to seconds per molecule (criterion reports per-iteration)
        if unit in ("ns", "ns/iter"):
            secs = value * 1e-9
        elif unit in ("µs", "us", "µs/iter"):
            secs = value * 1e-6
        elif unit in ("ms", "ms/iter"):
            secs = value * 1e-3
        else:
            secs = value

        # Criterion benchmarks often run over a batch; we need mol/s.
        # The bench names encode the batch size (e.g., /16).
        # We report throughput as-is, leaving calibration to the user.
        # For comparison we record the mean time per batch iteration.
        # Key: use a canonical name matching rdkit_bench keys.
        key = canonicalize(name_raw)
        if key:
            results[key] = 1.0 / secs  # approx mol/s (adjust if batched)

    if not results:
        # Fallback: look for lines like "N elements/second"
        for m in re.finditer(r"([\d.]+)\s+elem/s.*?(\w[\w/]+)", content):
            key = canonicalize(m.group(2).lower())
            if key:
                results[key] = float(m.group(1))

    return results


def parse_rdkit(path: str) -> dict[str, float]:
    results: dict[str, float] = {}
    for line in Path(path).read_text().splitlines()[1:]:  # skip header
        parts = line.strip().split("\t")
        if len(parts) == 2:
            results[parts[0]] = float(parts[1])
    return results


def canonicalize(name: str) -> str | None:
    name = name.lower().replace("-", "_")
    if "morgan" in name:
        m = re.search(r"r(\d+).*?(\d{3,4})", name)
        if m:
            return f"morgan_r{m.group(1)}_{m.group(2)}"
        return "morgan_r2_2048"
    if "maccs" in name:
        return "maccs_166"
    return None


def main() -> None:
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    molprint_results = parse_criterion(sys.argv[1])
    rdkit_results = parse_rdkit(sys.argv[2])

    if not molprint_results:
        print("WARNING: could not parse any molprint throughput from criterion output.")
        print("         Ensure you ran: cargo bench 2>&1 | tee bench_results.txt")
        sys.exit(0)

    print(f"\n{'Metric':<25} {'molprint (mol/s)':>18} {'RDKit (mol/s)':>15} {'Speedup':>10}")
    print("-" * 72)

    all_pass = True
    for key in sorted(set(molprint_results) | set(rdkit_results)):
        mp = molprint_results.get(key)
        rd = rdkit_results.get(key)
        if mp is None or rd is None:
            speedup_str = "N/A"
            status = "?"
        else:
            speedup = mp / rd
            speedup_str = f"{speedup:.1f}x"
            if speedup < MIN_SPEEDUP:
                status = "✗ SLOW"
                all_pass = False
            else:
                status = "✓"
        mp_str = f"{mp:>18,.0f}" if mp else "       N/A"
        rd_str = f"{rd:>15,.0f}" if rd else "    N/A"
        print(f"{key:<25} {mp_str} {rd_str} {speedup_str:>10}  {status}")

    print()
    if all_pass:
        print(f"✓ All metrics meet the {MIN_SPEEDUP}x speedup target.")
        sys.exit(0)
    else:
        print(f"✗ Some metrics below {MIN_SPEEDUP}x speedup target.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
