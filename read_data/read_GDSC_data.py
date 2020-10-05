# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2019/12/18

READ GDSC DATA

NOTES:
- 2020/01/29: Add read_one_GDSC_modality for incorporating
and concatenating more tissues.
"""

import os
import pandas as pd
import numpy as np
from functools import reduce


def read_GDSC_data(data_types,
                   tissues,
                   projects=None,
                   folder='../../data/GDSC/'):
    
    # If folder is None, set to default
    data_folder = folder or '../../data/GDSC/'
    
    if tissues is None:
        tissues = 'all'
    if type(tissues) == str:
        tissues = [tissues]
    if type(data_types) == str:
        data_types = [data_types]
        
    data_df = {'%s_%s'%('-'.join(tissues),d): read_one_GDSC_modality(tissues, projects, d, data_folder)
               for d in data_types}
    
    # Reduce to same samples
    common_samples = reduce(np.intersect1d, [list(data_df[d].index.get_level_values(0).astype(str)) for d in data_df])
    data_df = {d:data_df[d].loc[common_samples] for d in data_df}
    assert(np.unique([data_df[d].shape[0] for d in data_df]).shape[0] == 1)
    
    return data_df


def read_one_GDSC_modality(tissues, projects, data_type, folder):
    if type(tissues) != list:
        tissues = [tissues]
    if type(projects) != list:
        projects = [projects]
    
    folders = ['%s/%s/GDSC_%s_data.pkl'%(folder,
                                        data_type,
                                        data_type)
              for t in tissues for p in projects]
    folders = np.array(folders)
    
    folders = folders[[os.path.exists(f) for f in folders]]
    if folders.shape[0] == 0:
        return pd.DataFrame()

    data_df = [pd.read_pickle(f, compression='gzip') for f in folders]

    # Take same features
    unique_features = reduce(np.intersect1d, [df.columns for df in data_df])
    data_df = pd.concat([df[unique_features] for df in data_df])
    return data_df