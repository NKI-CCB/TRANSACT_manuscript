import os
import pandas as pd

import read_downloaded_files
import aggregate_gdc_data

# Transform data into design matrix
file_annot = read_downloaded_files.read_metadata('./TCGA_manifests/manifest.txt')
data_df, features_df = aggregate_gdc_data.aggregate_gdc_data('rnaseq', '../data/TCGA/raw/')

# Biospecimen
biospecimen_file = '../data/TCGA/biospecimen.txt'
biospecimen_df = pd.read_csv(biospecimen_file, sep=',')

# Merge read counts and biospecimen to get barcode
file_annot = file_annot.merge(biospecimen_df, how='left', left_on='aliquot_id', right_on='aliquot')
data_df = data_df.merge(file_annot[['id', 'barcode']].set_index('id'),
                     how='left',
                     left_index=True,
                     right_index=True)
data_df = data_df.set_index('barcode')

# Save data as GZIP pickle
data_df.to_pickle('../data/TCGA/rnaseq/TCGA_rnaseq_data.pkl', compression='gzip')
file_annot.to_pickle('../data//TCGA/rnaseq/TCGA_rnaseq_sample_annot.pkl', compression='gzip')
features_df.to_pickle('../data/TCGA/rnaseq/TCGA_rnaseq_feature_annot.pkl', compression='gzip')