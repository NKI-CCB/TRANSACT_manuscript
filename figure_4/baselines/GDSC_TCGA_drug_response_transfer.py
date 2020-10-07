#!/usr/bin/env python
# coding: utf-8

## Drug response study
#
# For a list of drugs:
#   - Align cell line data and tumor.
#   - Project cell line data and train predictor.
#   - Use predictor to predict on tumors.
#   - Compare prediction. 

import os, sys

######
#
# PARAMETERS
#
######

tissues = {
    'GDSC': ['All'],
    'TCGA': ['All']
}
projects = {
    'GDSC': [None],
    'TCGA': [None]
}

data_sources = ['GDSC', 'TCGA']

data_types = ['rnaseq']
genes_filtering = 'mini_cancer'
data_normalization = 'library_size' # Can be TPM, "library_size" or "log". Else will not have any influence.

source = 'GDSC'
target = 'TCGA'

with_mean = True
with_std = True

if len(sys.argv) > 1 and sys.argv[1] == 'linear':
    kernel_surname = 'PRECISE'
    kernel_name = 'linear'
    kernel_param = {}
elif len(sys.argv) > 1 and sys.argv[1] == 'rbf':
    kernel_surname = 'TRANSACT'
    kernel_name = 'rbf'
    kernel_param = {
        'gamma': 0.0005
    }
else:
    kernel_surname = 'default'
    kernel_name = 'rbf'
    kernel_param = {
        'gamma': 0.0005
    }

print(kernel_surname)

number_pc = {
    'source': 70,
    'target': 150
}

n_pv = 30
n_interpolation = 100

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
    ('Doxorubicin', 133, 'Doxorubicin'),
    ('Docetaxel', 1819, 'Docetaxel'),
    ('Cyclophosphamide', None, 'Cyclophosphamide'),
    ('Etoposide', None, 'Etoposide'),
    ('Oxaliplatin', 1806, 'Oxaliplatin'),
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


from transact.TRANSACT import TRANSACT


#######
# Folders for figures
#######

# General data folder
figure_folder = './figures/'
kernel_subfolder = kernel_surname
if kernel_subfolder in os.listdir(figure_folder):
    print('BEWARE: ALREADY COMPUTATION IN FIGURE FILE')
else:
    os.makedirs(figure_folder + kernel_subfolder)

kernel_subfolder = figure_folder + kernel_subfolder

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
# Alignment
#######


# Compute principal vectors
TRANSACT_clf = TRANSACT(kernel=kernel_name,
                        kernel_params=kernel_param,
                        n_components=number_pc,
                        n_jobs=20,
                        verbose=10)

TRANSACT_clf.fit(normalized_data_df[source_data_key],
                normalized_data_df[target_data_key],  
                n_pv=n_pv,
                step=n_interpolation,
                with_interpolation=True)

# Project data
source_consensus_features = TRANSACT_clf.transform(normalized_data_df[source_data_key])
target_consensus_features = TRANSACT_clf.transform(normalized_data_df[target_data_key])

# Put into format
source_consensus_features = pd.DataFrame(source_consensus_features,
                                         index=data_df[source_data_key].index)
target_consensus_features = pd.DataFrame(target_consensus_features,
                                         index=data_df[target_data_key].index)


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
    
    # Train predictor
    TRANSACT_clf.fit_predictor(X_source, y_source.values.flatten())

    # Predict value
    y_target_predicted = TRANSACT_clf.predict(X_target)
    y_target_df = pd.DataFrame.copy(y_target)
    y_target_df['predicted'] = np.array(y_target_predicted).astype(float)

    # Merge response data
    y_target_df = y_target_df.replace('Clinical Progressive Disease', 'Non Responder')
    y_target_df = y_target_df.replace('Stable Disease', 'Non Responder')
    y_target_df = y_target_df.replace('Partial Response', 'Responder')
    y_target_df = y_target_df.replace('Complete Response', 'Responder')
    order = ['Responder', 'Non Responder']
    box_pairs = [('Responder', 'Non Responder')]

    # Plot
    plt.figure(figsize=(5,7))
    bplot = sns.boxplot(x='measure_of_response',
                        y='predicted',
                        data=y_target_df,
                        order=order,
                        linewidth=2.,
                        width=.8,
                        whis=[5,95],
                        showfliers=False,
                        boxprops=dict(alpha=.4))

    sns.swarmplot(x='measure_of_response',
                y='predicted',
                data=y_target_df,
                order=order)

    annot = add_stat_annotation(bplot, data=y_target_df,
                                x='measure_of_response', y='predicted',
                                box_pairs = box_pairs,
                                order=order,
                                text_format='full',
                                test=test,
                                comparisons_correction=None,
                                loc='inside',
                                verbose=2,
                                fontsize=16)
    plt.xlabel('')
    plt.xticks(np.arange(len(order)), [e.replace(' ', '\n') for e in order], fontsize=20, color='black')
    plt.yticks(fontsize=15, color='black')
    plt.ylabel('Predicted AUC', fontsize=25, color='black')
    plt.tight_layout()

    file_save = '%s/predicted_AUC_n_pv_%s_GDSC_%s_TCGA_%s%s_%s.png'%(kernel_subfolder,
                                                                    n_pv,
                                                                    GDSC_drug_name,
                                                                    TCGA_drug_name,
                                                                    '_%s'%(GDSC_drug_id) if GDSC_drug_id else '',
                                                                    test)
    plt.savefig(file_save, dpi=300)
    pd.DataFrame(annot[1]).to_csv(file_save.replace('.png', '.csv'))
    with open(file_save.replace('.png', '_sample_size.txt'), 'w') as size_file :
        size_file.write('GDSC,%s,\nTCGA,%s,'%(X_source.shape[0], X_target.shape[0]))
    with open(file_save.replace('.png', '_effect_size.txt'), 'w') as es_file:
        responders_val = y_target_df[y_target_df['measure_of_response'] == 'Responder']['predicted'].values
        non_responders_val = y_target_df[y_target_df['measure_of_response'] == 'Non Responder']['predicted'].values
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