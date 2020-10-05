# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2020/06/17

READ PDXE DRUG RESPONSE DATA

"""

import os
import pandas as pd
import numpy as np
from functools import reduce

def read_PDXE_response(PDXE_drug_response_df, PDXE_drug_name, X_target):
    # X_target has to be a DataFrame with genes in columns

    y_target = PDXE_drug_response_df[PDXE_drug_response_df['Treatment'] == PDXE_drug_name]
    y_target = y_target.set_index('Model')
    
    if np.unique([len(x) for x in X_target.index])[0] == 2:
        X_target.index = [e[0] for e in X_target.index]

    common_samples = np.intersect1d(list(y_target.index), X_target.index)
    y_target = y_target.loc[common_samples]
    X_target_response = X_target.loc[common_samples]

    del common_samples
    return X_target_response, y_target