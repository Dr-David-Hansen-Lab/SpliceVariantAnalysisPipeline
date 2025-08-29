#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

# Simple and fast GTF parser for transcript->gene mapping
# Assumes attributes contain gene_id and transcript_id

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
    rows = []
    with open(gtf_path, 'r') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != 'transcript':
                # Some GTFs annotate transcript info on exon lines as well; collect from any line with transcript_id
                attrs = parse_attributes(parts[8])
                if 'transcript_id' in attrs and 'gene_id' in attrs:
                    rows.append((attrs['transcript_id'], attrs['gene_id'], attrs.get('gene_name', attrs['gene_id'])))
                continue
            attrs = parse_attributes(parts[8])
            if 'transcript_id' in attrs and 'gene_id' in attrs:
                rows.append((attrs['transcript_id'], attrs['gene_id'], attrs.get('gene_name', attrs['gene_id'])))
    # Deduplicate
    df = pd.DataFrame(rows, columns=['transcript_id', 'gene_id', 'gene_name']).drop_duplicates()
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
