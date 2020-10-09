#!/usr/bin/env python
# coding: utf-8

## ComBat correction between GDSC and PDX followed by GDSC cross-validation
#
# In this script, I reproduce the pipeline from [Sakelaropoulos et al 2019] 
# for GDSC-PDXE analysis. After correction, I consider one drug and compute 
# for various hyper-parameters of the neural networks the predictive 
# performance by 5-fold cross-validation.
# The process goes as follows:
#   - we compute a KFold to divide dataset in 5.
#   - for each fold held out and each parameter possible, we train on the
#   remaining folds and predict the response.
#   - the final measure (either pearson correlation or MSE) is computed among
#   the combined held-out samples.

n_splits = 5
output_folder = './output/GDSC_to_TCGA/'
figure_folder = './figures/GDSC_TCGA_combat/'

import os, sys, getopt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import functools
from joblib import Parallel, delayed
import scipy
from datetime import date
import re
import uuid
from pickle import dump
from combat.pycombat import pycombat
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

from clf_utils import make_network, make_figure_folder
from read_GDSC import read_GDSC_drug_response


# possible names: Erlotinib, Gemcitabine, Cetuximab, Afatinib, Paclitaxel, Alpelisib, 5-Fluorouracil
#.                Trametinib, Ruxolitinib, Ribociclib, Tamoxifen, Docetaxel, ...
GDSC_drug_name = None
GDSC_drug_id = None # Or 1007/1819 for Docetaxel
random_state = None

opts, args = getopt.getopt(sys.argv[1:],"n:i:r:",["name=","id="])
for opt, arg in opts:
    if opt in ("-i", "--ifile"):
        GDSC_drug_id = int(arg)
    elif opt in ("-n", "--ofile"):
        GDSC_drug_name = arg
    elif opt in ("-r", "--ofile"):
        random_state = int(arg)

random_state = random_state if random_state else 15485

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

with_mean = False
with_std = False


## Import data
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

# remove some genes to avoid ComBat to collapse
number_top_genes = 1700

top_source_variable_genes = pd.DataFrame(np.var(data_df[source_data_key]), columns=['variance'])
top_source_variable_genes = top_source_variable_genes.sort_values('variance', ascending=False)
top_source_variable_genes = top_source_variable_genes.head(number_top_genes).index
top_target_variable_genes = pd.DataFrame(np.var(data_df[target_data_key]), columns=['variance'])
top_target_variable_genes = top_target_variable_genes.sort_values('variance', ascending=False)
top_target_variable_genes = top_target_variable_genes.head(number_top_genes).index
top_variable_genes = np.intersect1d(top_source_variable_genes, top_target_variable_genes)

for d in data_df:
    data_df[d] = data_df[d][top_variable_genes]


## Correct with ComBat
data_corrected = pycombat(pd.concat(list(data_df.values())).T,
                         [1]*data_df[source_data_key].shape[0] + [2]*data_df[target_data_key].shape[0])

normalized_data_df = {
    k : data_corrected[data_df[k].index].T
    for k in data_df
}
normalized_data_df[source_data_key].index = pd.MultiIndex.from_tuples(normalized_data_df[source_data_key].index)


# Read response
if GDSC_drug_name in ['Cetuximab',
                    'Doxorubicin',
                    'Etoposide',
                    'Bleomycin',
                    'Bicalutamide',
                    'Bleomycin (50 uM)', 
                    'Pemetrexed',
                    'AICA Ribonucleotide']:
    GDSC_drug_response_file = '../data/GDSC/response/GDSC1_drug_response.xlsx'
else:
    GDSC_drug_response_file = '../data/GDSC/response/GDSC1_drug_response.xlsx'
GDSC_drug_response_df = pd.read_excel(GDSC_drug_response_file)

X, y = read_GDSC_drug_response(GDSC_drug_name,
                            GDSC_drug_id,
                            GDSC_drug_response_df,
                            normalized_data_df[source_data_key])

y_values = y.values.astype(np.float32)
X_values = X.values.astype(np.float32)

## Cross-validation
# Splitting data
split_folds = KFold(n_splits=5, shuffle=True, random_state=random_state)
split_folds.get_n_splits()


# Import parameters
param_dict = load(open('./params/sakellaropoulos_param_grid.pkl', 'rb'))
parameter_grid = ParameterGrid(param_dict)

if GDSC_drug_name not in os.listdir(output_folder):
    os.mkdir(output_folder + GDSC_drug_name)
output_folder += GDSC_drug_name + '/'

# Cross-validation for each parameter in the grid
for param_idx, param in enumerate(parameter_grid):
    print('START PARAM NUMBER %s'%(param_idx))

    param['n_input'] = X_values.shape[1]
    network_folder = make_figure_folder(output_folder, param)

    y_predict_nn = np.zeros(y_values.shape[0])

    for split_idx, (train_index, valid_index) in enumerate(split_folds.split(X_values, y_values)):
        print('START SPLIT NUMBER %s'%(split_idx))
        X_train, y_train = X_values[train_index], y_values[train_index]
        X_valid, y_valid = X_values[valid_index], y_values[valid_index]
        
        net = make_network(param)
        net = NeuralNetRegressor(
            net,
            max_epochs=param['n_epochs'],
            lr=param['learning_rate'],
            batch_size=param['batch_size'],
            device= 'cuda' if torch.cuda.is_available() else 'cpu',
            optimizer=torch.optim.SGD,
            optimizer__momentum=param['momentum'],
            optimizer__weight_decay=param['l2_penalty'],
            iterator_train__shuffle = True,
            verbose=0
        )
        pipeline = Pipeline([
            ('scaler', StandardScaler(with_mean=with_mean, with_std=with_std)),
            ('net', net)
        ]
        )
        pipeline.fit(X_train,
                     y_train)
        y_predict_nn[valid_index] = pipeline.predict(X_valid.astype(np.float32)).flatten()
        
        pipeline_file = 'clf_random-state_%s_split-nbre_%s.pkl'%(random_state,
                                                                split_idx)
        dump(pipeline, '%s/%s'%(network_folder, pipeline_file))

    pd.DataFrame(
        y_predict_nn,
        index=y.index
    ).to_csv('%s/prediction_random-state_%s.csv'%(network_folder,
                                                  random_state))

    pred_perf = scipy.stats.pearsonr(y_predict_nn, y_values.flatten())
    MSE = mean_squared_error(y_predict_nn, y_values.flatten())
    perf_df = pd.DataFrame([pred_perf[0], MSE], index=['pred_perf', 'MSE']).T
    perf_df.to_csv('%s/pred_perf_random-state_%s.csv'%(network_folder,
                                                       random_state))

    dump(param, '%s/param.pkl'%(network_folder))


# Cross-validation with ElasticNet
alpha_values = np.logspace(-5,10,16)
l1_ratio_values = np.linspace(0,10,11)/10
param_grid ={
    'regression__alpha': alpha_values,
    'regression__l1_ratio': l1_ratio_values
}

y_predict_baseline = np.zeros(y_values.shape[0])
for split_idx, (train_index, valid_index) in enumerate(split_folds.split(X_values, y_values)):
    X_train, y_train = X_values[train_index], y_values[train_index]
    X_valid, y_valid = X_values[valid_index], y_values[valid_index]
    
    baseline_grid = GridSearchCV(Pipeline([
                            ('scaler', StandardScaler(with_mean=with_mean, with_std=with_std)),
                            ('regression', ElasticNet())
                            ]),
                                 cv=10,
                                 n_jobs=30,
                                 param_grid=param_grid,
                                 verbose=1,
                                 scoring='neg_mean_squared_error')
    baseline_grid.fit(X_train, y_train)
    y_predict_baseline[valid_index] = baseline_grid.predict(X_valid)

pd.DataFrame(
        y_predict_baseline,
        index=y.index
).to_csv('%s/baseline_prediction_random-state_%s.csv'%(output_folder,
                                                    random_state))

pred_perf_baseline = scipy.stats.pearsonr(y_predict_baseline, y_values.flatten())
MSE_baseline = mean_squared_error(y_predict_baseline, y_values.flatten())
perf_df = pd.DataFrame([pred_perf_baseline[0], MSE_baseline], index=['pred_perf', 'MSE']).T
perf_df.to_csv('%s/baseline_pred_perf_random-state_%s.csv'%(output_folder,
                                                            random_state))


## Results processing
fig_name = '%s%s_%s%s_%s'%(figure_folder,
                           GDSC_drug_name,
                           '_centered' if with_mean else '',
                           '_standardized' if with_std else '',
                           random_state)

param_names = ['hidden', 'input', 'activation', 'hiddenDO', 'inputDO', 'l2pen', 'lr']
def parse_folder_results(f):
    param = {}
    for n in param_names:
        param[n] = re.search('%s_([0-9A-Za-z-.]+)'%(n), f)
        param[n] = [param[n].group(1)] if param[n] else ''
    param['folder'] = f
    param_df = pd.DataFrame.from_dict(param)
    
    results_files = ['%s/%s/'%(output_folder, f) + e for e in os.listdir('%s/%s'%(output_folder, f))
                     if '.csv' in e and 'pred_perf' in e and (str(random_state) in e or random_state is None)]
    
    if len(results_files) == 0:
        return None
    
    results_df = [pd.read_csv(r, header=0, index_col=0) for r in results_files]
    results_df = pd.concat(results_df)
    results_df.index = [f] * results_df.shape[0]
        
    return results_df


relevant_subfolders = [e for e in os.listdir(output_folder)
                       if 'hidden' in e]

results_df = [parse_folder_results(f) for f in relevant_subfolders]
results_df = [df for df in results_df if df is not None]
results_df = pd.concat(results_df)

baseline_df = pd.read_csv('%s/baseline_pred_perf_random-state_%s.csv'%(output_folder, random_state),
                          header=0, index_col=0)

results_df.columns = [('model', e) for e in results_df.columns]
for e in ['MSE', 'pred_perf']:
    results_df[('baseline', e)] = baseline_df[e].values[0]
results_df.columns = pd.MultiIndex.from_tuples(results_df.columns)

results_df.to_csv('%s.csv'%(fig_name))

for x in ['pred_perf', 'MSE']:
    plt.figure(figsize=(10,4.5))
    plt.axhline(baseline_df[x].values[0],
                linewidth=3,
                label='ElasticNet',
                linestyle='dashed',
                color='black',
                alpha=0.6)
    plt.plot(results_df['model'].sort_values(x, ascending=(x=='pred_perf'))[x].values, '+',
             label='Neural network', markersize=10)
    plt.xticks(fontsize=15, color='black')
    plt.xlabel('Cross-validated models', fontsize=20, color='black')
    plt.yticks(fontsize=15, color='black')
    plt.ylabel('Predictive performance' if x == 'pred_perf' else 'Mean Squared Error', fontsize=20, color='black')
    plt.legend(fontsize=15, bbox_to_anchor=(1., 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('%s_CV_%s.png'%(fig_name, x),
                dpi=300)
    plt.clf()
    
print('%s models computed'%(results_df.shape[0]))