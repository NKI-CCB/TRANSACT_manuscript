# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2019/12/18

READ GDSC DRUG RESPONSE DATA

"""
import os
import pandas as pd
import numpy as np
from functools import reduce


def read_GDSC_response(GDSC_drug_response_frames, GDSC_drug_name, X_source, GDSC_drug_id=None):
    # X_source has to be a DataFrame with genes in columns

    if GDSC_drug_name in ['Cetuximab',
                          'Doxorubicin',
                          'Etoposide',
                          'Bleomycin',
                          'Bicalutamide',
                          'Bleomycin (50 uM)', 
                          'Pemetrexed',
                          'AICA Ribonucleotide']:
        GDSC_drug_response_df = GDSC_drug_response_frames['GDSC1'].copy()
    else:
        GDSC_drug_response_df = GDSC_drug_response_frames['GDSC2'].copy()
    
    y_source = GDSC_drug_response_df[GDSC_drug_response_df['DRUG_NAME'] == GDSC_drug_name]
    if GDSC_drug_id is not None:
        y_source = GDSC_drug_response_df[GDSC_drug_response_df['DRUG_ID'] == GDSC_drug_id]
    else:
        print(np.unique(GDSC_drug_response_df['DRUG_ID']).shape)
    y_source = y_source[['CELL_LINE_NAME', 'AUC']]
    y_source = y_source.set_index('CELL_LINE_NAME')

    while type(X_source.index) == pd.core.indexes.multi.MultiIndex:
        X_source.index = X_source.index.droplevel(1)
    X_source_response = X_source.groupby(X_source.index).mean()

    common_samples = np.intersect1d(y_source.index,
                                    X_source_response.index)
    y_source = y_source.loc[common_samples]
    X_source_response = X_source_response.loc[common_samples]

    return X_source_response, y_source