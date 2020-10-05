# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2020/01/23

READ DATA

Data reader from various sources with correction for same features.
This script outputs a dictionary with dataset from different sources and
different genomic views. Each view has the same features in the exact 
same order across the different data sources.
"""

import sys
import pandas as pd
import numpy as np
from functools import reduce

from read_GDSC_data import read_GDSC_data
from read_TCGA_data import read_TCGA_data
from read_PDXE_data import read_PDXE_data
from read_HMF_data import read_HMF_data


upload_methods = {
    'GDSC': read_GDSC_data,
    'TCGA': read_TCGA_data,
    'PDXE': read_PDXE_data,
    'HMF': read_HMF_data
}

def read_data(tissues,
              data_types,
              projects=None,
              data_sources=['GDSC', 'TCGA'],
              folder_basis='../../data'):
    """
    Read data for one type

    Parameters
    ----------
    tissues : string or list
        Tissues to consider. ex: ['Breast']

    data_types : string or list
        Data-types to consider. Implemented: 'fpkm', 'rnaseq', 'mut', 'cnv'.

    projects: string or list
        Projects to use, ex: ['BRCA']

    data_source: list, default to ['GDSC', 'BRCA']
        Name of data sources. Implemented: 'PDXE', 'FPKM', 'GDSC'.

    folder: dictionary, default to None
        Name of the folders where the data must be located

    Returns
    -------
    data_df: dictionary
        Dictionary with data.
        Keys as data-source, values as dictionary {data-type: pd.DataFrame} with samples in rows.
        Rows are harmonized for each sub-dictionary (same sample order).
        Columns have same order across data-sources for the same data-type.
    """
    
    # Read data
    if type(tissues) != dict:
        tissues = {
            s:tissues for s in data_sources
        }
    if type(projects) != dict:
        projects = {
            s:projects for s in data_sources
        }

    folders = {}
    for s in data_sources:
        folders[s] = folder_basis + s

    data_df = {
        s: upload_methods[s](data_types=data_types,
                             tissues=tissues[s],
                             projects=projects[s],
                             folder=folders[s])
        for s in data_sources
    }
    print('DATA UPLOADED')
                
    # Harmonize to same features
    unique_data_types = np.unique([list(data_df[s].keys()) for s in data_sources])
    for d in data_types:
        relevant_data_types = {s:[e for e in unique_data_types if d in e and e in data_df[s]]\
                               for s in data_sources}
        features = [data_df[s][e].columns for s in data_sources for e in relevant_data_types[s]]
        unique_features = reduce(np.intersect1d, features)
        for s in data_sources:
            for e in relevant_data_types[s]:
                data_df[s][e] = data_df[s][e][unique_features]
                data_df[s][e] = data_df[s][e].groupby(by=data_df[s][e].columns, axis=1).agg('mean')
    print('GENE HARMONIZED ACROSS VIEWS')
    
    # Check good harmonization
    unique_data_identifiers = reduce(np.intersect1d, [list(data_df[s].keys()) for s in data_sources])
    for identif in unique_data_identifiers:
        features = [data_df[s][identif].columns for s in data_df]
        [np.testing.assert_array_equal(f,g) for f in features for g in features]
    print('CHECKED FOR GOOD HARMONIZATION')
        
    
    return data_df