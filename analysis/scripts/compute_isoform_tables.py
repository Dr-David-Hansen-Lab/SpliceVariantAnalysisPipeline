#!/usr/bin/env python3
import argparse
import os
import pandas as pd

"""
Reads a Salmon quant.sf and a tx2gene map and outputs:
- proportions table: transcript, gene, TPM, proportion
- dominant isoform table: gene, transcript, proportion
"""

def read_tx2gene(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t')
    # normalize columns
    cols = {c.lower(): c for c in df.columns}
    if 'transcript_id' not in cols and 'transcript' in cols:
        df.rename(columns={cols['transcript']: 'transcript_id'}, inplace=True)
    if 'gene_id' not in cols and 'gene' in cols:
        df.rename(columns={cols['gene']: 'gene_id'}, inplace=True)
    # Ensure expected names
    df.rename(columns={cols.get('transcript_id', 'transcript_id'): 'transcript_id',
                       cols.get('gene_id', 'gene_id'): 'gene_id'}, inplace=True)
    # Keep gene_name if present
    keep = ['transcript_id', 'gene_id'] + ([ 'gene_name' ] if 'gene_name' in df.columns else [])
    return df[keep]

def load_quant(quant_path: str) -> pd.DataFrame:
    q = pd.read_csv(quant_path, sep='\t')
    # Expect column Name and TPM; NumReads is preferred for counts plots if available
    if 'Name' not in q.columns or 'TPM' not in q.columns:
        raise ValueError('quant.sf missing required columns Name and TPM')
    # Salmon often encodes additional metadata in Name separated by '|'.
    # Extract the transcript_id as the first field before the first '|'.
    names = q['Name'].astype(str)
    # Fast split: if any '|' present, split, else keep as-is
    if names.str.contains('|', regex=False).any():
        transcript_id = names.str.split('|', n=1, expand=True)[0]
    else:
        transcript_id = names
    # Assemble output with TPM and NumReads if present
    cols = ['TPM'] + (['NumReads'] if 'NumReads' in q.columns else [])
    out = q[cols].copy()
    out.insert(0, 'transcript_id', transcript_id)
    return out

def compute_tables(quant: pd.DataFrame, tx2gene: pd.DataFrame):
    m = quant.merge(tx2gene, on='transcript_id', how='left')
    # handle transcripts not in mapping: set gene_id to NA and drop for proportions
    m_valid = m.dropna(subset=['gene_id']).copy()
    # avoid divide-by-zero: if gene sum TPM == 0, set proportion to 0
    sums = m_valid.groupby('gene_id')['TPM'].transform('sum')
    m_valid['proportion'] = 0.0
    nonzero = sums > 0
    m_valid.loc[nonzero, 'proportion'] = m_valid.loc[nonzero, 'TPM'] / sums[nonzero]
    cols = ['transcript_id', 'gene_id', 'TPM', 'proportion']
    # Carry NumReads through if present
    if 'NumReads' in m_valid.columns:
        cols.insert(3, 'NumReads')
    if 'gene_name' in m_valid.columns:
        cols.insert(2, 'gene_name')
    props = m_valid[cols]
    # dominant isoform per gene by highest proportion
    idx = props.groupby('gene_id')['proportion'].idxmax()
    dcols = ['gene_id', 'transcript_id', 'proportion']
    if 'gene_name' in props.columns:
        dcols.insert(1, 'gene_name')
    dom = props.loc[idx][dcols].reset_index(drop=True)
    return props, dom

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quant', required=True, help='Path to Salmon quant.sf')
    ap.add_argument('--tx2gene', required=True, help='Path to tx2gene TSV')
    ap.add_argument('--props', required=True, help='Output CSV for proportions')
    ap.add_argument('--dominant', required=True, help='Output CSV for dominant isoforms')
    args = ap.parse_args()

    tx2gene = read_tx2gene(args.tx2gene)
    quant = load_quant(args.quant)
    props, dom = compute_tables(quant, tx2gene)

    os.makedirs(os.path.dirname(args.props), exist_ok=True)
    props.to_csv(args.props, index=False)
    dom.to_csv(args.dominant, index=False)

if __name__ == '__main__':
    main()
