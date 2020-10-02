import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce


# FPKM
fpkm_file = '../data/PDXE/fpkm/PDXE_fpkm_data.pkl'
fpkm_df = pd.read_pickle(fpkm_file, compression='gzip')
fpkm_df.columns = [e.split('.')[0] for e in fpkm_df.columns]

# Mini cancer
mini_cancer_df = pd.read_csv('../data/mini_cancer_genes.csv', index_col=0)

# Save
common_genes = np.intersect1d(fpkm_df.columns, mini_cancer_df['Hugo'])
fpkm_df = fpkm_df[common_genes]

fpkm_df.to_pickle(fpkm_file, compression='gzip')