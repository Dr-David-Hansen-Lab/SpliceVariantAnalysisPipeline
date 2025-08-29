#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

"""
Reads all sample proportion CSVs from a directory and plots per-gene stacked barplots
of isoform proportions, one PNG per gene.
"""

def load_props(props_dir: str) -> pd.DataFrame:
    dfs = []
    for fn in os.listdir(props_dir):
        if not fn.endswith('_proportions.csv'):
            continue
        sample = fn.replace('_proportions.csv', '')
        df = pd.read_csv(os.path.join(props_dir, fn))
        df['sample'] = sample
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError('No *_proportions.csv files found in ' + props_dir)
    return pd.concat(dfs, ignore_index=True)


def plot_gene(df: pd.DataFrame, gene: str, out_dir: str, metric: str = 'proportion'):
    # Accept either gene_id or gene_name matches
    if 'gene_name' in df.columns:
        gdf = df[(df['gene_id'] == gene) | (df['gene_name'] == gene)].copy()
    else:
        gdf = df[df['gene_id'] == gene].copy()
    if gdf.empty:
        # create empty placeholder figure
        plt.figure(figsize=(6, 3))
        plt.text(0.5, 0.5, f'No data for {gene}', ha='center', va='center')
        plt.axis('off')
        os.makedirs(out_dir, exist_ok=True)
        plt.savefig(os.path.join(out_dir, f'{gene}_proportions.png'), dpi=200, bbox_inches='tight')
        plt.close()
        return

    # Make a stacked bar chart: x=sample, hue=transcript_id, weight=metric (proportion or NumReads)
    # Pivot to wide and plot
    value_col = metric if metric in gdf.columns else 'proportion'
    pivot = gdf.pivot_table(index='sample', columns='transcript_id', values=value_col, aggfunc='sum', fill_value=0)
    pivot = pivot.sort_index()

    os.makedirs(out_dir, exist_ok=True)
    plt.figure(figsize=(max(6, len(pivot) * 0.8), 4))
    import numpy as np
    bottom = None
    colors = sns.color_palette('tab20', n_colors=max(3, pivot.shape[1]))
    for i, col in enumerate(pivot.columns):
        vals = pivot[col].to_numpy(dtype=float)
        x = np.arange(len(pivot.index))
        if bottom is None:
            plt.bar(x, vals, label=col, color=colors[i])
            bottom = vals.copy()
        else:
            plt.bar(x, vals, bottom=bottom, label=col, color=colors[i])
            bottom = bottom + vals
    plt.ylabel('Isoform ' + ('proportion' if value_col == 'proportion' else 'counts'))
    plt.xlabel('Sample')
    title_gene = gene
    if 'gene_name' in gdf.columns and not gdf.empty:
        # Prefer showing gene_name (first non-null) with ID in parentheses if different
        nm = str(gdf['gene_name'].dropna().unique()[0]) if gdf['gene_name'].notna().any() else gene
        gid = str(gdf['gene_id'].dropna().unique()[0]) if gdf['gene_id'].notna().any() else gene
        title_gene = nm if nm == gid else f"{nm} ({gid})"
    title_metric = 'proportions' if value_col == 'proportion' else 'counts'
    plt.title(f'Isoform {title_metric} for {title_gene}')
    plt.xticks(ticks=np.arange(len(pivot.index)), labels=[str(x) for x in pivot.index.tolist()], rotation=45, ha='right')
    plt.legend(title='Transcript', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{gene}_proportions.png'), dpi=200)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--props_dir', required=True)
    ap.add_argument('--genes', nargs='+', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--metric', choices=['proportion', 'NumReads'], default='proportion', help='Value to plot: proportion or NumReads')
    args = ap.parse_args()

    df = load_props(args.props_dir)
    for gene in args.genes:
        plot_gene(df, gene, args.out_dir, metric=args.metric)

if __name__ == '__main__':
    main()
