import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce


# Read counts
rnaseq_file = '../data/TCGA/rnaseq/TCGA_rnaseq_data.pkl'
rnaseq_df = pd.read_pickle(rnaseq_file, compression='gzip')
rnaseq_df.columns = [e.split('.')[0] for e in rnaseq_df.columns]

# Mini cancer
mini_cancer_df = pd.read_csv('../data/mini_cancer_lookup_genes.csv', index_col=0)
mini_cancer_df = mini_cancer_df.reset_index().set_index('ENSEMBL')

# Save
common_genes = np.intersect1d(rnaseq_df.columns, mini_cancer_df.index)
common_mini_genes = mini_cancer_df.loc[common_genes]['Hugo'].values
assert common_genes.shape[0] == common_genes.shape[0]
rnaseq_df = rnaseq_df[common_genes]
rnaseq_df.columns = common_mini_genes

rnaseq_df.to_pickle(rnaseq_file, compression='gzip')

# Response data
response_file = '../data/TCGA/raw/bioinfo16_supplementary_tables.xlsx'
response_df = pd.read_excel(response_file, sheet_name=None)
response_df = response_df['Table S2']
response_df.columns = response_df.loc[1]
response_df = response_df.drop([0,1,2])
response_df.to_csv('../data/TCGA/response/response.csv')
