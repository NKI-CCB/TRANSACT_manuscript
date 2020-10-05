# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2020/07/05

READ TCGA DRUG RESPONSE DATA

"""
import os
import pandas as pd
import numpy as np
from functools import reduce

def read_TCGA_response(TCGA_drug_name, X_target, TCGA_drug_response_file):
    # Read data
    TCGA_drug_response_df = pd.read_csv(TCGA_drug_response_file, sep=',', header=[0,1])
    TCGA_drug_response_df.columns = TCGA_drug_response_df.columns.droplevel(1)

    # Restrict to one drug
    TCGA_drug_response_df = TCGA_drug_response_df[TCGA_drug_response_df['drug_name'] == TCGA_drug_name]
    TCGA_drug_response_df = TCGA_drug_response_df[['bcr_patient_barcode', 'drug_name', 'measure_of_response']]
    TCGA_drug_response_df = TCGA_drug_response_df.drop_duplicates()

    # Shorten barcodes
    length_barcode_input = np.min([len(e) for e in X_target.index])
    length_barcode_response = np.min([len(e) for e in TCGA_drug_response_df['bcr_patient_barcode'].values])
    length_barcode = min(length_barcode_input, length_barcode_response)

    TCGA_drug_response_df['bcr_patient_barcode'] = [e[:length_barcode] for e in TCGA_drug_response_df['bcr_patient_barcode']]
    y_target = TCGA_drug_response_df[['bcr_patient_barcode', 'measure_of_response']].set_index('bcr_patient_barcode')
    unique_patients_response = np.unique(y_target.index, return_counts=True)
    unique_patients_response = unique_patients_response[0][unique_patients_response[1] == 1]
    y_target = y_target.loc[unique_patients_response]
    X_target_response = X_target.copy()
    X_target_response.index = [e[:length_barcode] for e in X_target_response.index]

    # Take common samples
    common_TCGA_samples = np.intersect1d(X_target_response.index, y_target.index)
    y_target = y_target.loc[common_TCGA_samples]
    X_target_response = X_target_response.loc[common_TCGA_samples]
    X_target_response = X_target_response.groupby(X_target_response.index).agg('mean') #For redundant samples

    assert X_target_response.shape[0] == y_target.shape[0]
    assert np.all(X_target_response.index == y_target.index)

    return X_target_response, y_target