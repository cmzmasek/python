#!/usr/bin/env python3
"""
Parse all FASTA files in a directory and write the first N sequences to an output file.

Usage:
    python subset_fasta.py -d <input_dir> -n <num_sequences> -o <output_file>
    python subset_fasta.py -d . -n 10 -o first10.fasta
"""

import argparse
import sys
from pathlib import Path


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            elif line:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


HELP = """\
subset_fasta.py  —  extract the first N sequences from FASTA files in a directory
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

USAGE
    python subset_fasta.py -n <N> -o <output> [options]

REQUIRED
    -n, --num   <int>    Number of sequences to write to the output file
    -o, --output <file>  Output FASTA file path

OPTIONAL
    -d, --dir    <dir>    Input directory to search (default: current directory)
    -s, --suffix <ext>    File suffix to match (default: .fasta)
                          The leading dot is optional: fasta and .fasta both work
    -h, --help            Show this help message and exit

NOTES
    Files are processed in alphabetical order. Sequences are collected across
    files until N total have been written. No external dependencies required.

EXAMPLES
    # First 10 sequences from FASTA files in the current directory
    python subset_fasta.py -n 10 -o first10.fasta

    # First 50 sequences from a specific directory
    python subset_fasta.py -d /data/genomes -n 50 -o subset.fasta

    # Match files ending in .fa instead of .fasta
    python subset_fasta.py -n 25 -o out.fasta -s .fa

    # Leading dot is optional
    python subset_fasta.py -n 25 -o out.fasta -s fa
"""


def main():
    if len(sys.argv) == 1:
        print(HELP)
        sys.exit(0)

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("-d", "--dir", default=".")
    parser.add_argument("-n", "--num", type=int)
    parser.add_argument("-o", "--output")
    parser.add_argument("-s", "--suffix", default=".fasta")
    args = parser.parse_args()

    if args.help:
        print(HELP)
        sys.exit()

    errors = []
    if args.num is None:
        errors.append("  -n / --num is required")
    if args.output is None:
        errors.append("  -o / --output is required")
    if errors:
        print("Error: missing required argument(s):")
        print("\n".join(errors))
        print("\nRun with -h for help.")
        sys.exit(1)

    suffix = args.suffix if args.suffix.startswith(".") else f".{args.suffix}"

    input_dir = Path(args.dir).resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"Error: '{input_dir}' is not a directory.")

    fasta_files = sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and p.suffix.lower() == suffix.lower()
    )

    if not fasta_files:
        raise SystemExit(f"No FASTA files found in '{input_dir}' with suffix '{suffix}'")

    written = 0
    with open(args.output, "w") as out:
        for fasta_file in fasta_files:
            count = 0
            for header, seq in parse_fasta(fasta_file):
                if count >= args.num:
                    break
                out.write(f"{header}\n{seq}\n")
                count += 1
            written += count

    print(f"Wrote {written} sequence(s) from {len(fasta_files)} file(s) [{suffix}] to '{args.output}'.")


if __name__ == "__main__":
    main()
