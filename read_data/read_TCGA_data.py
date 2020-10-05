# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2020/01/16

READ TCGA DATA
"""

import os
import pandas as pd
import numpy as np
from functools import reduce


def read_one_TCGA_modality(tissues,
                    projects,
                    data_type,
                    folder):
    
    if type(tissues) != list:
        tissues = [tissues]
    if type(projects) != list:
        projects = [projects]

    folders = np.array(['%s/%s/TCGA_%s_data.pkl'%(folder,
                                                data_type,
                                                data_type)
              for t in tissues for p in projects])

    folders = folders[[os.path.exists(f) for f in folders]]
    if folders.shape[0] == 0:
        return pd.DataFrame()

    data_df = [pd.read_pickle(f, compression='gzip') for f in folders]

    # Take same features
    unique_features = reduce(np.intersect1d, [df.columns for df in data_df])
    data_df = pd.concat([df[unique_features] for df in data_df])
    return data_df


def read_TCGA_data(data_types=['rnaseq'],
                tissues=None, 
                projects=None,
                folder='/DATA/s.mourragui/data/2020_01_TCGA_data'):

    # If folder is None, set to default
    data_folder = folder or '/DATA/s.mourragui/data/2020_01_TCGA_data'

    if type(tissues) == str:
        tissues = [tissues]

    data_df = {'%s_%s'%('-'.join(tissues),d):read_one_TCGA_modality(tissues, projects, d, data_folder)\
            for d in data_types}
    data_df = {d:data_df[d] for d in data_df if not data_df[d].empty}
    
    for i in data_df:
        data_df[i].index = [e[:19] for e in data_df[i].index]
        # Merge samples since correlations are very good
        data_df[i] = data_df[i].reset_index().groupby('index').agg('mean') 
    unique_samples = reduce(np.intersect1d, [data_df[i].index for i in data_df])
    
    return {d:data_df[d].loc[unique_samples] for d in data_df}