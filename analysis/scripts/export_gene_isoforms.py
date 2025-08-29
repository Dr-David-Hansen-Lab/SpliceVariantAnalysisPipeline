#!/usr/bin/env python3
import argparse
import os
import pandas as pd

"""
Collect per-sample isoform tables and export gene-specific combined CSVs with TPM, NumReads, and proportion
for the provided genes (match on gene_id or gene_name).
Output columns: sample, gene_id, gene_name, transcript_id, TPM, NumReads, proportion
"""

def load_props_dir(props_dir: str) -> pd.DataFrame:
    dfs = []
    for fn in os.listdir(props_dir):
        if not fn.endswith('_proportions.csv'):
            continue
        sample = fn.replace('_proportions.csv', '')
        df = pd.read_csv(os.path.join(props_dir, fn))
        df['sample'] = sample
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError(f'No *_proportions.csv found in {props_dir}')
    return pd.concat(dfs, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--props_dir', required=True)
    ap.add_argument('--genes', nargs='+', required=True)
    ap.add_argument('--out_dir', required=True)
    args = ap.parse_args()

    df = load_props_dir(args.props_dir)
    os.makedirs(args.out_dir, exist_ok=True)

    for gene in args.genes:
        if 'gene_name' in df.columns:
            gdf = df[(df['gene_id'] == gene) | (df['gene_name'] == gene)].copy()
        else:
            gdf = df[df['gene_id'] == gene].copy()
        if gdf.empty:
            # write an empty file with headers for consistency
            cols = ['sample', 'gene_id', 'gene_name', 'transcript_id', 'TPM', 'NumReads', 'proportion']
            pd.DataFrame(columns=cols).to_csv(os.path.join(args.out_dir, f'{gene}_isoforms.csv'), index=False)
            continue
        # Ensure all expected columns exist
        for c in ['gene_name', 'NumReads']:
            if c not in gdf.columns:
                gdf[c] = pd.NA
        out_cols = ['sample', 'gene_id', 'gene_name', 'transcript_id', 'TPM', 'NumReads', 'proportion']
        gdf[out_cols].sort_values(['sample', 'NumReads'], ascending=[True, False]).to_csv(
            os.path.join(args.out_dir, f'{gene}_isoforms.csv'), index=False
        )


if __name__ == '__main__':
    main()
