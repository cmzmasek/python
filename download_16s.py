#!/usr/bin/env python3
"""
Download one 16S rRNA sequence per bacterial species from GenBank.

Input file: one species name per line, e.g.:
    Lactobacillus acidophilus
    Escherichia coli
    Bacillus subtilis

Usage:
    python download_16s.py species.txt --out sequences.fasta --email you@example.com
"""

import argparse
import os
import sys
import textwrap
import time
from Bio import Entrez, SeqIO

VERSION = '1.0.0'


def fetch_16s(species, email, n=1, api_key=None):
    """Search GenBank for up to n full-length 16S sequences for the given species."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    query = (
        f'"{species}"[Organism] AND '
        '(16S ribosomal RNA[Title] OR 16S rRNA[Title] OR 16S rDNA[Title]) AND '
        '1000:2000[SLEN]'
    )

    handle = Entrez.esearch(db='nucleotide', term=query, retmax=n, sort='relevance')
    record = Entrez.read(handle)
    handle.close()

    if not record['IdList']:
        return []

    ids = ','.join(record['IdList'])
    handle = Entrez.efetch(db='nucleotide', id=ids, rettype='fasta', retmode='text')
    seq_records = list(SeqIO.parse(handle, 'fasta'))
    handle.close()

    return seq_records


BANNER = f"""
  ┌─────────────────────────────────────────────────────────────┐
  │            16S rRNA Sequence Downloader  (GenBank)          │
  │                        version {VERSION:<28}│
  └─────────────────────────────────────────────────────────────┘

  Reads a list of bacterial species names and downloads one
  full-length 16S rRNA sequence per species from NCBI GenBank.
  Results are saved as a multi-FASTA file. Species with no hit
  are logged separately.
"""

EPILOG = textwrap.dedent("""\
  input file format:
    One species name per line. Blank lines and lines starting
    with '#' are ignored.

      Lactobacillus acidophilus
      Escherichia coli
      # this line is a comment
      Bacillus subtilis

  examples:
    python download_16s.py species.txt
    python download_16s.py -t "Escherichia coli" "Bacillus subtilis"
    python download_16s.py species.txt -t "Clostridioides difficile"
    python download_16s.py species.txt --email you@example.com
    python download_16s.py species.txt -n 3 --out results.fasta
    python download_16s.py species.txt --out results.fasta --missing no_hits.txt
    python download_16s.py species.txt --api-key ABC123

  notes:
    - Sequences are filtered to 1000–2000 bp (full-length 16S).
    - NCBI requires a valid email address for Entrez queries.
    - Get a free NCBI API key at: https://www.ncbi.nlm.nih.gov/account/
      (raises rate limit from 3 to 10 requests/second)
""")


def main():
    parser = argparse.ArgumentParser(
        prog='download_16s.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=BANNER,
        epilog=EPILOG,
        add_help=False,
    )

    parser.add_argument('input', nargs='?',
                        help='Text file with one species name per line')
    parser.add_argument('-t', dest='taxa', nargs='+', metavar='TAXON',
                        help='One or more species names to query directly')
    parser.add_argument('--out', default='16s_sequences.fasta', metavar='FILE',
                        help='Output FASTA file  (default: 16s_sequences.fasta)')
    parser.add_argument('--email', default='user@example.com', metavar='ADDR',
                        help='Your email address, required by NCBI  (default: user@example.com)')
    parser.add_argument('--api-key', dest='api_key', metavar='KEY',
                        help='NCBI API key for 10 req/s instead of 3')
    parser.add_argument('-n', type=int, default=1, metavar='N',
                        help='Max sequences to download per species  (default: 1)')
    parser.add_argument('--missing', default='missing.txt', metavar='FILE',
                        help='File to record species with no hits  (default: missing.txt)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Allow overwriting existing output files')
    parser.add_argument('-h', '--help', action='store_true',
                        help='Show this help message and exit')

    args = parser.parse_args()

    if args.help or (not args.input and not args.taxa):
        parser.print_help()
        sys.exit(0)

    if not args.overwrite:
        conflicts = [f for f in (args.out, args.missing) if os.path.exists(f)]
        if conflicts:
            for f in conflicts:
                print(f'Error: file already exists: {f}', file=sys.stderr)
            print('Use --overwrite to allow overwriting.', file=sys.stderr)
            sys.exit(1)

    species_list = []
    if args.input:
        with open(args.input) as f:
            species_list += [line.strip() for line in f if line.strip() and not line.startswith('#')]
    if args.taxa:
        species_list += args.taxa

    print(f'Processing {len(species_list)} species...', file=sys.stderr)

    found = []
    missing = []
    delay = 0.11 if args.api_key else 0.34  # stay within NCBI rate limits

    for i, species in enumerate(species_list, 1):
        print(f'[{i}/{len(species_list)}] {species}', file=sys.stderr, end=' ')
        try:
            records = fetch_16s(species, args.email, n=args.n, api_key=args.api_key)
            if records:
                for record in records:
                    record.description = f'[query: {species}] {record.description}'
                found.extend(records)
                ids = ', '.join(r.id for r in records)
                print(f'-> {ids}', file=sys.stderr)
            else:
                missing.append(species)
                print('-> NOT FOUND', file=sys.stderr)
        except Exception as e:
            print(f'-> ERROR: {e}', file=sys.stderr)
            missing.append(species)

        time.sleep(delay)

    with open(args.out, 'w') as f:
        SeqIO.write(found, f, 'fasta')

    if missing:
        with open(args.missing, 'w') as f:
            f.write('\n'.join(missing) + '\n')
        print(f'\nMissing sequences ({len(missing)}) written to {args.missing}', file=sys.stderr)

    print(f'\nDone: {len(found)} sequences saved to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
