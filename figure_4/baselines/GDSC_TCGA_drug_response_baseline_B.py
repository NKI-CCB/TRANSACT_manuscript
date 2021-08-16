#!/usr/bin/env python
# coding: utf-8

## DRUG RESPONSE ON TCGA FOR BASELINE B (Deep Learning)
#
# For a list of drugs:
#   - Look-up best set of hyper-parameter (trained on cell lines)
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
    ('Gemcitabine', None, 'Gemcitabine'),
    ('Paclitaxel', None, 'Paclitaxel'),
    ('Oxaliplatin', 1806, 'Oxaliplatin'),
    ('Afatinib', None, 'Trastuzumab'),
    ('Vinorelbine', None, 'Vinorelbine'),
    ('5-Fluorouracil', None, 'Fluorouracil'),
    ('Temozolomide', None, 'Temozolomide'),
    ('Doxorubicin', 133, 'Doxorubicin'),
    ('Docetaxel', 1819, 'Docetaxel'),
    ('Cyclophosphamide', None, 'Cyclophosphamide'),
    ('Etoposide', None, 'Etoposide'),
    ('Bleomycin', None, 'Bleomycin'),
    ('Pemetrexed', None, 'Pemetrexed'),
    ('Irinotecan', None, 'Irinotecan'),
    ('Cetuximab', None, 'Cetuximab')
]

output_uncorrected_cv_folder = '../data/output/GDSC_only/'
random_state = 183627362
n_retrain = 50

# Import 

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy
import uuid
from pickle import dump
from datetime import date
sns.set_style("whitegrid")
sns.set_context('paper')

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, KFold, GroupKFold, GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import ParameterGrid
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error
from sklearn.pipeline import Pipeline
from sklearn.utils import shuffle, resample
from joblib import dump, load, Parallel, delayed
from statannot.statannot import add_stat_annotation
import torch
from torch.autograd import Variable
import torch.nn.functional as F
import torch.utils.data as Data
from torch.utils.data import Dataset, TensorDataset, DataLoader
from torch.utils.data.dataset import random_split
from skorch import NeuralNetClassifier, NeuralNetRegressor

sys.path.insert(0, '../read_data/')
from read_data import read_data
from read_GDSC_response import read_GDSC_response
from read_PDXE_response import read_PDXE_response
from read_TCGA_response import read_TCGA_response
from reformat_df import reformat_df
import library_size_normalization

sys.path.insert(0, '../src/')
from clf_utils import make_network, make_drug_folder
from cv_results_processing import parse_folder_results, read_best_param, make_skorch_network

#####
#
# DATA FOLDER
#
#####
figure_folder = './output/'
figure_subfolder = 'baseline_B'
if figure_subfolder not in os.listdir(figure_folder):
    os.makedirs(figure_folder + figure_subfolder)
figure_folder += figure_subfolder + '/'

#####
#
# IMPORT DATA
#
#####
data_df = read_data(tissues=tissues,
                    data_types=[e for e in data_types],
                    projects=projects,
                    data_sources=data_sources,
                    folder_basis='../data/')

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

# Removing the healthy samples
healthy_samples_index = data_df[target_data_key].index.str.contains(r'-(10A|11A)-')
data_df[target_data_key] = data_df[target_data_key].loc[~healthy_samples_index]

for ds in list(data_df.keys()):
    GE_normalized = library_size_normalization.TMM_normalization(data_df[ds].values.astype(float))
    GE_normalized = np.array(GE_normalized)
    GE_normalized = np.log(np.array(GE_normalized)+1)
    data_df[ds] = pd.DataFrame(GE_normalized,
                               columns=data_df[ds].columns,
                               index=data_df[ds].index)

normalized_data_df = {
    k: StandardScaler(with_mean=with_mean, with_std=with_std).fit_transform(data_df[k])
    for k in data_df
}
normalized_data_df = {
    k: pd.DataFrame(normalized_data_df[k], index=data_df[k].index, columns=data_df[k].columns)
    for k in data_df
}

#####
#
# LOAD DRUG RESPONSE DATA
#
#####
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

for GDSC_drug_name, GDSC_drug_id, TCGA_drug_name in n_retrain * drug_list:
    print(GDSC_drug_name, TCGA_drug_name)
    drug_folder = make_drug_folder(GDSC_drug_name, TCGA_drug_name, figure_folder)
    training_id = str(uuid.uuid4())[:8]

    # Read best param by cross-val
    uncorrected_cv_folder = output_uncorrected_cv_folder + GDSC_drug_name +\
                            ('_centered' if with_mean else '') +\
                            ('_standardized' if with_std else '')
    uncorrected_param = read_best_param(uncorrected_cv_folder, random_state, '%s/%s/uncorrected_cv_results.csv'%(figure_folder, drug_folder))
    uncorrected_param['n_input'] = normalized_data_df[source_data_key].shape[1]
    
    # Import data
    X_source, y_source = read_GDSC_response(GDSC_drug_response_frames, 
                                        GDSC_drug_name,
                                        normalized_data_df[source_data_key].copy(),
                                        GDSC_drug_id)
    X_target, y_target = read_TCGA_response(TCGA_drug_name,
                                            normalized_data_df[target_data_key].copy(),
                                            TCGA_drug_response_file)

    # Train network
    net = make_network(uncorrected_param)
    net = NeuralNetRegressor(
        net,
        max_epochs=uncorrected_param['n_epochs'],
        lr=uncorrected_param['learning_rate'],
        batch_size=uncorrected_param['batch_size'],
        device= 'cuda' if torch.cuda.is_available() else 'cpu',
        optimizer=torch.optim.SGD,
        train_split=None,
        optimizer__momentum=uncorrected_param['momentum'],
        optimizer__weight_decay=uncorrected_param['l2_penalty'],
        iterator_train__shuffle = True,
        verbose=0
    )
    pipeline = Pipeline([
        ('net', net)
    ])
    pipeline.fit(X_source.values.astype(np.float32),
                 y_source.values.astype(np.float32))
    dump(pipeline, open('%s/%s/clf_%s.pkl'%(figure_folder, drug_folder, training_id), 'wb'))

    # Predict on TCGA
    y_target_predicted_raw = pipeline.predict(X_target.values.astype(np.float32))
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

    file_save = '%s/%s/baseline_predicted_AUC_ID_%s%s_%s.png'%(figure_folder,
                                                            drug_folder,
                                                            training_id,
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