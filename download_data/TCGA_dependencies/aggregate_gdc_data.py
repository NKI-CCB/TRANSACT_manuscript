# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui
with code from Tycho Bismeijer

2020/01/03

AGGREGATE GDC DATA

Aggregate GDC downloaded data
"""

import os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from read_downloaded_files import read_one_methyl_file, read_one_rnaseq_file, read_one_miRNA_file

def aggregate_gdc_data(data_type, raw_folder):
    if data_type == 'methyl':
        reading_method = read_one_methyl_file
    elif data_type in ['rnaseq', 'fpkm']:
        reading_method = read_one_rnaseq_file
    elif data_type.lower() == 'mirna':
        reading_method = read_one_miRNA_file

    data_df = Parallel(n_jobs=1, verbose=1)(delayed(reading_method)(s, raw_folder)
                                                  for s in os.listdir(raw_folder))

    # Remove files that do not have same length that the other ones
    data_df_shapes = [r[1].shape[1] for r in data_df if r]
    data_df_shapes = np.unique(data_df_shapes, return_counts=True)
    consensual_shape = data_df_shapes[0][np.argmax(data_df_shapes[1])]
    data_df = [r for r in data_df if r and r[1].shape[1] == consensual_shape]

    # Test that same features are being used.
    if data_type == 'methyl':
        columns_of_interest = ['Chromosome', 'Start', 'End']
        for i,(_,_,f) in enumerate(data_df):
            pd.testing.assert_frame_equal(data_df[0][-1][columns_of_interest],f[columns_of_interest])
    if data_type == 'miRNA':
        columns_of_interest = ['reads_per_million_miRNA_mapped', 'miRNA_ID']
        data_df = [(r[0], r[1][columns_of_interest].set_index('miRNA_ID').T, r[2]) for r in data_df]

    features_df = data_df[0][-1]

    # Concatenate data
    data_df = list(zip(*[(a,b) for a,b,_ in data_df]))
    data_index = data_df[0]
    if data_type in ['rnaseq', 'fpkm']:
        data_df = pd.concat([x.set_index('gene').T for x in data_df[1]])
    elif data_type in ['methyl', 'miRNA']:
        data_df = pd.concat([x for x in data_df[1]])
        
    data_df.index = data_index
    
    return data_df, features_df