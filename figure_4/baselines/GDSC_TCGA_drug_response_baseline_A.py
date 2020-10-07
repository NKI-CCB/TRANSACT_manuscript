#!/usr/bin/env python
# coding: utf-8

## Drug response baseline
#
# For a list of drugs:
#   - Train on cell lines.
#   - Apply on tumors.
# This is a simple baseline to compare the approach.

######
#
# PARAMETERS
#
######

tissues = {
    'GDSC': ['All'],
    'TCGA': ['TCGA']
}
projects = {
    'GDSC': [None],
    'TCGA': ['all']
}

data_sources = ['GDSC', 'TCGA']

data_types = ['rnaseq']
genes_filtering = 'mini_cancer'
data_normalization = 'library_size' # Can be TPM, "library_size" or "log". Else will not have any influence.

source = 'GDSC'
target = 'TCGA'

with_mean = True
with_std = True

test = 'Mann-Whitney-ls'

# Order: (GDSC_drug_name, GDSC_drug_ID, TCGA_drug_name)
drug_list =[
    ('Cisplatin', None, 'Cisplatin'),
    ('Cisplatin', None, 'Carboplatin'),
    ('Afatinib', None, 'Trastuzumab'),
    ('Gemcitabine', None, 'Gemcitabine'),
    ('Paclitaxel', None, 'Paclitaxel'),
    ('5-Fluorouracil', None, 'Fluorouracil'),
    ('Temozolomide', None, 'Temozolomide'),
    ('Doxorubicin', 1386, 'Doxorubicin'),
    ('Docetaxel', 1819, 'Docetaxel'),
    ('Cyclophosphamide', None, 'Cyclophosphamide'),
    ('Etoposide', None, 'Etoposide'),
    ('Oxaliplatin', 1806, 'Oxaliplatin'),
    ('Bleomycin (50 uM)', None, 'Bleomycin'),
    ('Bleomycin', None, 'Bleomycin'),
    ('Pemetrexed', None, 'Pemetrexed'),
    ('Vinorelbine', None, 'Vinorelbine'),
    ('Irinotecan', None, 'Irinotecan'),
    ('Cetuximab', None, 'Cetuximab')
]


# Import 

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy
from datetime import date
sns.set_style("whitegrid")
sns.set_context('paper')

from sklearn.model_selection import KFold, GridSearchCV
from sklearn.linear_model import Ridge, ElasticNet, Lasso
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from statannot.statannot import add_stat_annotation

sys.path.insert(0, '../read_data/')
from read_data import read_data
from read_GDSC_response import read_GDSC_response
from read_PDXE_response import read_PDXE_response
from read_TCGA_response import read_TCGA_response
from reformat_df import reformat_df
import library_size_normalization


#######
# Folders for figures
#######

# General data folder
figure_folder = './figures/'
figure_subfolder = 'baseline_A'
if figure_subfolder not in os.listdir(figure_folder):
    os.makedirs(figure_folder + figure_subfolder)
figure_folder += figure_subfolder

#######
# Import data
#######

# Import and format data

data_df = read_data(tissues=tissues,
                    data_types=[e for e in data_types],
                    projects=projects,
                    data_sources=data_sources,
                    folder_basis='../data/')

source_data_key, target_data_key = reformat_df(data_df, source, target)

# Removing the healthy samples
healthy_samples_index = data_df[target_data_key].index.str.contains(r'-(10A|11A)-')
data_df[target_data_key] = data_df[target_data_key].loc[~healthy_samples_index]

# Library size correction
average_depth_global = 10**5
for ds in list(data_df.keys()):
    GE_normalized = library_size_normalization.TMM_normalization(data_df[ds].values.astype(float))
    GE_normalized = np.array(GE_normalized)
    average_depths = np.mean(np.sum(GE_normalized,1))
    GE_normalized = GE_normalized / average_depths * average_depth_global
    GE_normalized = np.log(np.array(GE_normalized)+1)
    data_df[ds] = pd.DataFrame(GE_normalized,
                               columns=data_df[ds].columns,
                               index=data_df[ds].index)


# Normalize data
normalized_data_df = {
    ds : StandardScaler(with_mean=with_mean, with_std=with_std).fit_transform(data_df[ds])
    for ds in data_df
}
for ds in normalized_data_df:
    normalized_data_df[ds] = pd.DataFrame(normalized_data_df[ds],
                                         index=data_df[ds].index,
                                         columns=data_df[ds].columns)


#######
# Drug response
#######:

# Load data
# GDSC
unique_drugs = None
GDSC_drug_response_frames = {}
for x in ['GDSC2', 'GDSC1']:
    GDSC_drug_response_file = '../data/GDSC/response/%s_fitted_dose_response_25Feb20.xlsx'%(x)
    GDSC_drug_response_frames[x] = pd.read_excel(GDSC_drug_response_file)
    if unique_drugs is None:
        unique_drugs = np.unique(GDSC_drug_response_frames[x]['DRUG_NAME'])
    else:
        unique_drugs = np.concatenate([unique_drugs, np.unique(GDSC_drug_response_frames[x]['DRUG_NAME'])])

# TCGA
TCGA_drug_response_file = '../data/TCGA/response/response.csv'


def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)



for GDSC_drug_name, GDSC_drug_id, TCGA_drug_name in drug_list:
    X_source, y_source = read_GDSC_response(GDSC_drug_response_frames, 
                                        GDSC_drug_name,
                                        normalized_data_df[source_data_key].copy(),
                                        GDSC_drug_id)
    X_target, y_target = read_TCGA_response(TCGA_drug_name,
                                            normalized_data_df[target_data_key].copy(),
                                            TCGA_drug_response_file)
    
    # Fit Elastic Net model
    param_grid ={
        'regression__alpha': np.logspace(-5,15,26)
    }
    grid_raw = GridSearchCV(Pipeline([('regression', ElasticNet(0.5))]),
                            cv=10, 
                            n_jobs=10, 
                            param_grid=param_grid, 
                            verbose=1, 
                            scoring='neg_mean_squared_error')
    grid_raw.fit(X_source, y_source)

    y_target_predicted_raw = grid_raw.predict(X_target)
    y_target_df = pd.DataFrame.copy(y_target)
    y_target_df['predicted raw'] = np.array(y_target_predicted_raw).astype(float)

    # Merge response data
    y_target_df = y_target_df.replace('Clinical Progressive Disease', 'Non Responder')
    y_target_df = y_target_df.replace('Stable Disease', 'Non Responder')
    y_target_df = y_target_df.replace('Partial Response', 'Responder')
    y_target_df = y_target_df.replace('Complete Response', 'Responder')
    order = ['Responder', 'Non Responder']
    box_pairs = [('Responder', 'Non Responder')]

    # Plot and save
    plt.figure(figsize=(5,7))
    bplot = sns.boxplot(x='measure_of_response',
                        y='predicted raw',
                        data=y_target_df,
                        order=order,
                        linewidth=2.,
                        width=.8,
                        whis=[5,95],
                        showfliers=False,
                        boxprops=dict(alpha=.4))

    sns.swarmplot(x='measure_of_response',
                y='predicted raw',
                data=y_target_df,
                order=order)

    annot_raw = add_stat_annotation(bplot, data=y_target_df,
                                    x='measure_of_response', y='predicted raw',
                                    box_pairs=box_pairs,
                                    order=order,
                                    text_format='full',
                                    test=test,
                                    loc='inside',
                                    verbose=2,
                                    fontsize=16)

    plt.xlabel('')
    plt.xticks(np.arange(len(order)), [e.replace(' ', '\n') for e in order], fontsize=20, color='black')
    plt.yticks(fontsize=15, color='black')
    plt.ylabel('Predicted AUC', fontsize=25, color='black')
    plt.tight_layout()

    file_save = '%s/baseline_predicted_AUC_GDSC_%s_TCGA_%s%s_%s.png'%(figure_folder,
                                                                        GDSC_drug_name,
                                                                        TCGA_drug_name,
                                                                        '_%s'%(GDSC_drug_id) if GDSC_drug_id else '',
                                                                       test)
    plt.savefig(file_save, dpi=300)
    pd.DataFrame(annot_raw[1]).to_csv(file_save.replace('.png', '.csv'))
    with open(file_save.replace('.png', '_sample_size.txt'), 'w') as size_file :
        size_file.write('GDSC,%s,\nTCGA,%s,'%(X_source.shape[0], X_target.shape[0]))
    with open(file_save.replace('.png', '_effect_size.txt'), 'w') as es_file:
        responders_val = y_target_df[y_target_df['measure_of_response'] == 'Responder']['predicted raw'].values
        non_responders_val = y_target_df[y_target_df['measure_of_response'] == 'Non Responder']['predicted raw'].values
        cohen_stat = cohen_d(non_responders_val, responders_val)
        total_effect = np.mean(non_responders_val) - np.mean(responders_val)
        es_file.write('cohen,%s\nmean_diff,%s'%(cohen_stat, total_effect))
    with open(file_save.replace('.png', '_ustat.txt'), 'w') as ustat_file:
        u_stat = scipy.stats.mannwhitneyu(non_responders_val, responders_val, alternative='greater')
        ustat_file.write('ustat,%s\npval,%s\nproduct_samples,%s'%(u_stat[0],
                                                                  u_stat[1],
                                                                  responders_val.shape[0]*non_responders_val.shape[0]))
        
    plt.clf()

    del file_save


