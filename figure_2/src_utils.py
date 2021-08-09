import pandas as pd
import re
import os
from skorch import NeuralNetClassifier, NeuralNetRegressor
import torch
import torch.optim

param_names = ['hidden', 'input', 'activation', 'hiddenDO', 'inputDO', 'l2pen', 'lr']

def parse_folder_results(f, folder, random_state):
    param = {}
    for n in param_names:
        param[n] = re.search('%s_([0-9A-Za-z-.]+)'%(n), f)
        param[n] = [param[n].group(1)] if param[n] else ''
    param['folder'] = f
    param_df = pd.DataFrame.from_dict(param)
    
    results_files = ['%s/%s/'%(folder, f) + e for e in os.listdir('%s/%s'%(folder, f))
                     if '.csv' in e and 'pred_perf' in e and (str(random_state) in e or random_state is None)]
    
    if len(results_files) == 0:
        return None
    
    results_df = [pd.read_csv(r, header=0, index_col=0) for r in results_files]
    results_df = pd.concat(results_df)
    results_df.index = [f] * results_df.shape[0]
        
    return results_df


def read_best_param(folder, random_state, output_fig=None, drug_folder_name=None):
    relevant_subfolders = [e for e in os.listdir(folder)
                           if 'hidden' in e]

    results_df = [parse_folder_results(f, folder, random_state)
                              for f in relevant_subfolders]
    results_df = [df for df in results_df if df is not None]
    results_df = pd.concat(results_df)

    baseline_df = pd.read_csv('%s/baseline_pred_perf_random-state_%s.csv'%(folder,
                                                                           random_state),
                              header=0, index_col=0)

    results_df.columns = [('model', e) for e in results_df.columns]
    for e in ['MSE', 'pred_perf']:
        results_df[('baseline', e)] = baseline_df[e].values[0]
    results_df.columns = pd.MultiIndex.from_tuples(results_df.columns)

    if output_fig is not None:
        results_df.to_csv('%s/%s'%(drug_folder_name, output_fig))
    
    best_model = results_df.sort_values(('model', 'pred_perf'), ascending=False).index[0]
    best_model_param = folder + '/' + best_model + '/param.pkl'
    best_model_param = load(open(best_model_param, 'rb'))
    return best_model_param


def make_skorch_network(net, param):
    return NeuralNetRegressor(
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