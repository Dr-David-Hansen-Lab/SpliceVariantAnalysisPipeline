"""
isoform_selection.py
Module for isoform selection and proportion calculation from Salmon quantification results.
"""
import pandas as pd
from typing import List, Dict

def load_quant_sf(quant_path: str) -> pd.DataFrame:
    """
    Load a Salmon quant.sf file into a DataFrame.
    Args:
        quant_path (str): Path to quant.sf file.
    Returns:
        pd.DataFrame: DataFrame with quantification data.
    """
    return pd.read_csv(quant_path, sep='\t')

def compute_isoform_proportions(quant_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute isoform proportions per gene from quant.sf DataFrame.
    Args:
        quant_df (pd.DataFrame): DataFrame with quant.sf data (must include 'Name', 'TPM', 'GeneID').
    Returns:
        pd.DataFrame: DataFrame with isoform proportions per gene.
    """
    # Group by gene and compute proportions
    quant_df['proportion'] = quant_df.groupby('GeneID')['TPM'].transform(lambda x: x / x.sum())
    return quant_df[['Name', 'GeneID', 'TPM', 'proportion']]

def get_dominant_isoform(quant_df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify the dominant isoform (highest proportion) per gene.
    Args:
        quant_df (pd.DataFrame): DataFrame with isoform proportions (must include 'GeneID', 'proportion').
    Returns:
        pd.DataFrame: DataFrame with dominant isoform per gene.
    """
    idx = quant_df.groupby('GeneID')['proportion'].idxmax()
    return quant_df.loc[idx][['GeneID', 'Name', 'proportion']]
