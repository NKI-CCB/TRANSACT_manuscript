import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce

# Read models
model_file = '../data/GDSC/model_list_20191104.csv'
model_df = pd.read_csv(model_file)

# Read counts
rnaseq_file = '../data/GDSC/raw/rnaseq_read_count_20191101.csv'
rnaseq_df = pd.read_csv(rnaseq_file, index_col=[0,1], header=[0,1,2])
rnaseq_df.columns = rnaseq_df.columns.droplevel(0).droplevel(1)
rnaseq_df.index = rnaseq_df.index.droplevel(0)
rnaseq_df = rnaseq_df.T
rnaseq_df = rnaseq_df.merge(model_df[['model_name', 'tissue']].set_index('model_name'),
                            how='left', left_index=True, right_index=True)
rnaseq_df = rnaseq_df.reset_index().set_index(['model_name', 'tissue'])
rnaseq_df = rnaseq_df.dropna(axis=1)

# FPKM
fpkm_file = '../data/GDSC/raw/rnaseq_fpkm_20191101.csv'

## NEW VERSION (2020/02)
fpkm_df = pd.read_csv(fpkm_file, index_col=[0,1], header=[0,1,2])
fpkm_df.columns = fpkm_df.columns.droplevel(0).droplevel(1)
fpkm_df.index = fpkm_df.index.droplevel(0)
fpkm_df = fpkm_df.T
fpkm_df = fpkm_df.dropna(axis=1)
fpkm_df = fpkm_df.merge(model_df[['model_name', 'tissue']].set_index('model_name'),
                        how='left', left_index=True, right_index=True)
fpkm_df = fpkm_df.reset_index().set_index(['model_name', 'tissue'])
fpkm_df = fpkm_df.dropna(axis=1)

# Mini cancer
mini_cancer_df = pd.read_csv('../data/mini_cancer_genes.csv', index_col=0)

# Save
common_genes_read = np.intersect1d(rnaseq_df.columns, mini_cancer_df['Hugo'])
rnaseq_df = rnaseq_df[common_genes_read]
rnaseq_df.to_pickle('../data/GDSC/rnaseq/GDSC_rnaseq_data.pkl', compression='gzip')

common_genes_fpkm = np.intersect1d(fpkm_df.columns, mini_cancer_df['Hugo'])
fpkm_df = fpkm_df[common_genes_read]
fpkm_df.to_pickle('../data/GDSC/fpkm/GDSC_fpkm_data.pkl', compression='gzip')

