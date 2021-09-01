#!/usr/bin/env python
# coding: utf-8

## Cross-validation of neural network on GDSC
#
# In this script, I consider one drug and compute for various hyper-parameters
# of the neural networks the predictive performance by 5-fold cross-validation.
# The process goes as follows:
#   - we compute a KFold to divide dataset in 5.
#   - for each fold held out and each parameter possible, we train on the
#   remaining folds and predict the response.
#   - the final measure (either pearson correlation or MSE) is computed among
#   the combined held-out samples.

with_mean = True
with_std = True
n_splits = 5
output_folder = './output/GDSC/'

import sys, getopt, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import functools
from joblib import Parallel, delayed
import scipy
import uuid
from pickle import load
sns.set_style("whitegrid")
sns.set_context('paper')

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, KFold, GroupKFold, GridSearchCV
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

from read_GDSC import read_GDSC, read_GDSC_drug_response
from clf_utils import split_data, torch_dataset_reduced, make_figure_folder, make_network, train_batch
from plot_perf import plot_pred_perf, plot_loss


device = 'cuda' if torch.cuda.is_available() else 'cpu'


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

# Read data
data_df = read_data(tissues={'GDSC':['All']},
                    data_types=['rnaseq'],
                    projects={'GDSC':[None]},
                    data_sources=['GDSC'],
                    folder_basis='../data/')

while type(data_df) == dict:
    keys = list(data_df.keys())
    if len(keys) != 1:
        print('WARNING: more than one data type')
        assert False
    data_df = data_df[keys[0]]

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
    GDSC_drug_response_file = '../data/GDSC/response/GDSC2_drug_response.xlsx'
GDSC_drug_response_df = pd.read_excel(GDSC_drug_response_file)

X, y = read_GDSC_drug_response(GDSC_drug_name,
                            GDSC_drug_id,
                            GDSC_drug_response_df,
                            data_df)

X_values = X.values
y_values = y.values.flatten()

drug_folder = '%s%s%s%s'%(GDSC_drug_name,
                        '_%s'%(GDSC_drug_id) if GDSC_drug_id else '',
                        '_centered' if with_mean else '',
                        '_standardized' if with_std else '')
if drug_folder not in os.listdir(output_folder):
    os.mkdir(output_folder + drug_folder)
output_folder += drug_folder + '/'


# Split data
split_folds = KFold(n_splits=5, shuffle=True, random_state=random_state)
split_folds.get_n_splits()


# Compute predictive performance for ElasticNet (baseline)
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

pred_perf_baseline = scipy.stats.pearsonr(y_predict_baseline, y_values)
MSE_baseline = mean_squared_error(y_predict_baseline, y_values)
perf_df = pd.DataFrame([pred_perf_baseline[0], MSE_baseline], index=['pred_perf', 'MSE']).T
perf_df.to_csv('%s/baseline_pred_perf_random-state_%s.csv'%(output_folder,
                                                            random_state))


# Import neural networks possible parameters
# Import parameters
param_dict = load(open('./params/sakellaropoulos_param_grid.pkl', 'rb'))
parameter_grid = ParameterGrid(param_dict)


for param_idx, param in enumerate(parameter_grid):
    print('START PARAM NUMBER %s'%(param_idx))

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
            iterator_train__shuffle = True
        )
        pipeline = Pipeline([
            ('scaler', StandardScaler(with_mean=with_mean, with_std=with_std)),
            ('net', net)
        ]
        )
        pipeline.fit(X_train.astype(np.float32),
                     y_train.reshape(-1,1).astype(np.float32))
        y_predict_nn[valid_index] = pipeline.predict(X_valid.astype(np.float32)).flatten()
        
        pipeline_file = 'clf_random-state_%s_split-nbre_%s.pkl'%(random_state,
                                                                split_idx)
        dump(pipeline, '%s/%s'%(network_folder, pipeline_file))

    pd.DataFrame(
        y_predict_nn,
        index=y.index
    ).to_csv('%s/prediction_random-state_%s.csv'%(network_folder,
                                                  random_state))

    pred_perf = scipy.stats.pearsonr(y_predict_nn, y_values)
    MSE = mean_squared_error(y_predict_nn, y_values)
    perf_df = pd.DataFrame([pred_perf[0], MSE], index=['pred_perf', 'MSE']).T
    perf_df.to_csv('%s/pred_perf_random-state_%s.csv'%(network_folder,
                                                       random_state))

    dump(param, '%s/param.pkl'%(network_folder))

