import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce

PDXE_folder = '../data/PDXE/'

raw_data_file = os.listdir('%sraw/'%(PDXE_folder))
raw_data_file = [f for f in raw_data_file if '.xlsx' in f]
assert len(raw_data_file) == 1
raw_data_file = raw_data_file[0]
raw_df = pd.read_excel('%s/raw/%s'%(PDXE_folder, raw_data_file), sheet_name=None)

fpkm_df = raw_df['RNAseq_fpkm']
pct_raw_df = raw_df['PCT raw data']
pct_curve_df = raw_df['PCT curve metrics']

# Save annotations
tissues_df = pct_raw_df[['Model', 'Tumor Type']]
tissues_df = tissues_df.drop_duplicates().dropna()
tissues_df = tissues_df.set_index('Model')
tissues_df.to_csv('%spancancer_biospecimen.csv'%(PDXE_folder), sep=',')

# Save FPKM
fpkm_df = fpkm_df.set_index('Sample').T
fpkm_df = fpkm_df.merge(tissues_df, left_index=True, right_index=True, how='left')
fpkm_df = fpkm_df.fillna('NA')
fpkm_df.to_pickle('%s/fpkm/PDXE_fpkm_data.pkl'%(PDXE_folder), compression='gzip')

# Save response
pct_curve_df.to_csv('../data/PDXE/response/response.csv')