# isoform_analysis.py
"""
Python module for extracting and summarizing isoform quantification results for user-specified genes.
Reads Salmon quant.sf files and outputs isoform proportions for selected genes.
"""
import pandas as pd
import os
from typing import List, Dict

def load_quant_sf(quant_path: str) -> pd.DataFrame:
    """
    Load a Salmon quant.sf file into a DataFrame.
    Args:
        quant_path: Path to quant.sf file.
    Returns:
        DataFrame with columns: Name, Length, EffectiveLength, TPM, NumReads
    """
    return pd.read_csv(quant_path, sep='\t')

def filter_isoforms_by_gene(quant_df: pd.DataFrame, gene_map: Dict[str, str], genes_of_interest: List[str]) -> pd.DataFrame:
    """
    Filter isoforms belonging to genes of interest.
    Args:
        quant_df: DataFrame from quant.sf
        gene_map: Dict mapping transcript IDs to gene IDs
        genes_of_interest: List of gene IDs or symbols
    Returns:
        Filtered DataFrame with isoforms for selected genes
    """
    quant_df['gene_id'] = quant_df['Name'].map(gene_map)
    return quant_df[quant_df['gene_id'].isin(genes_of_interest)]

def compute_isoform_proportions(filtered_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute isoform proportions per gene.
    Args:
        filtered_df: DataFrame with isoforms for selected genes
    Returns:
        DataFrame with isoform proportions per gene
    """
    result = filtered_df.copy()
    result['gene_total_TPM'] = result.groupby('gene_id')['TPM'].transform('sum')
    result['isoform_proportion'] = result['TPM'] / result['gene_total_TPM']
    return result

# Add more functions as needed for summarization, plotting, etc.
