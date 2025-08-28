"""
test_isoform_selection.py
Unit tests for isoform_selection.py
"""
import pandas as pd
import pytest
from analysis.isoform_selection import load_quant_sf, compute_isoform_proportions, get_dominant_isoform

def test_compute_isoform_proportions():
    data = {
        'Name': ['tx1', 'tx2', 'tx3', 'tx4'],
        'GeneID': ['geneA', 'geneA', 'geneB', 'geneB'],
        'TPM': [80, 20, 10, 90]
    }
    df = pd.DataFrame(data)
    result = compute_isoform_proportions(df)
    assert pytest.approx(result[result['Name']=='tx1']['proportion'].values[0], 0.01) == 0.8
    assert pytest.approx(result[result['Name']=='tx2']['proportion'].values[0], 0.01) == 0.2
    assert pytest.approx(result[result['Name']=='tx3']['proportion'].values[0], 0.01) == 0.1
    assert pytest.approx(result[result['Name']=='tx4']['proportion'].values[0], 0.01) == 0.9

def test_get_dominant_isoform():
    data = {
        'Name': ['tx1', 'tx2', 'tx3', 'tx4'],
        'GeneID': ['geneA', 'geneA', 'geneB', 'geneB'],
        'TPM': [80, 20, 10, 90],
        'proportion': [0.8, 0.2, 0.1, 0.9]
    }
    df = pd.DataFrame(data)
    dom = get_dominant_isoform(df)
    assert set(dom['Name']) == {'tx1', 'tx4'}
    assert set(dom['GeneID']) == {'geneA', 'geneB'}
