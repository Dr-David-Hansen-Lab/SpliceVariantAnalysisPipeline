#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

# Simple and fast GTF parser for transcript->gene mapping
# Assumes attributes contain gene_id and transcript_id
# Deduplicates on the fly using a set/dict to minimize memory and I/O

def parse_attributes(attr: str) -> dict:
    out = {}
    for field in attr.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        if ' ' in field:
            key, val = field.split(' ', 1)
            out[key] = val.strip().strip('"')
    return out

def build_tx2gene(gtf_path: str) -> pd.DataFrame:
    # Map (transcript_id, gene_id) -> gene_name (prefer explicit gene_name over fallback gene_id)
    best_name = {}
    with open(gtf_path, 'r') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            attrs = parse_attributes(parts[8])
            if 'transcript_id' not in attrs or 'gene_id' not in attrs:
                continue
            tx = attrs['transcript_id']
            gid = attrs['gene_id']
            gname = attrs.get('gene_name', gid)
            key = (tx, gid)
            # Insert once; upgrade to non-default gene_name if encountered later
            if key not in best_name:
                best_name[key] = gname
            else:
                if best_name[key] == gid and gname != gid:
                    best_name[key] = gname
    if not best_name:
        return pd.DataFrame(columns=['transcript_id', 'gene_id', 'gene_name'])
    rows = [(tx, gid, best_name[(tx, gid)]) for (tx, gid) in best_name.keys()]
    # Sort for reproducibility
    rows.sort(key=lambda x: (x[1], x[0]))  # by gene_id then transcript_id
    df = pd.DataFrame(rows, columns=['transcript_id', 'gene_id', 'gene_name'])
    return df


def main():
    p = argparse.ArgumentParser(description='Create tx2gene mapping from GTF')
    p.add_argument('--gtf', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()

    df = build_tx2gene(args.gtf)
    if df.empty:
        print('Warning: No transcript->gene mappings found in GTF', file=sys.stderr)
    df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()
