# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2020/01/16

REFORMAT DATA DF
"""

import pandas as pd
import numpy as np


def reformat_df(data_df, source, target):
    for ds in list(data_df.keys()):
        assert len(data_df[ds].keys()) == 1
        new_key = ('%s_%s'%(ds, list(data_df[ds].keys())[0])).replace('fpkm', 'tpm')
        data_df[new_key] = data_df[ds][list(data_df[ds].keys())[0]]
        del data_df[ds]
        
    source_data_key = [ds for ds in data_df if source in ds]
    assert len(source_data_key) == 1
    source_data_key = np.unique(source_data_key)[0]

    target_data_key = [ds for ds in data_df if target in ds]
    assert len(target_data_key) == 1
    target_data_key = np.unique(target_data_key)[0]

    return source_data_key, target_data_key