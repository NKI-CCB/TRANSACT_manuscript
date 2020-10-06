#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import uuid
import functools
import torch
import torch.utils.data as Data
from torch.utils.data import Dataset, TensorDataset, DataLoader
from torch.utils.data.dataset import random_split
from torch.autograd import Variable
import torch.nn.functional as F
import torch.utils.data as Data
from torch.utils.data import Dataset, TensorDataset, DataLoader
from torch.utils.data.dataset import random_split



def split_data(X_source, y_source, output_folder, proportion_train=0.8):
    X = torch.from_numpy(X_source.values).float()
    y = torch.from_numpy(y_source.values).float()
    data = TensorDataset(X, y)

    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, shuffle=True)
    train_size = int(proportion_train * X.shape[0])
    train_dataset, val_dataset = random_split(data, [train_size, X.shape[0]-train_size])
    train_loader = DataLoader(dataset=train_dataset, batch_size=20)
    val_loader = DataLoader(dataset=val_dataset, batch_size=X.shape[0]-train_size)

    split_uuid = str(uuid.uuid4())[:4]
    if 'split_idx' not in os.listdir(output_folder):
        os.mkdir('%s/split_idx'%(output_folder))
    np.savetxt('%s/split_idx/val_idx_%s.csv'%(output_folder, split_uuid),
                val_dataset.indices)

    return train_dataset, val_dataset, train_loader, val_loader, split_uuid


def torch_dataset_reduced(ds):
    # Reduce Torch dataset to numpy array by restricting indices
    X, y = ds.dataset.tensors
    return X[ds.indices], y[ds.indices]


def make_figure_folder(output_folder, param):
    # Make figure subfolder
    figure_subfolder = 'hidden_%s_input_%s_activation_%s_hiddenDO_%s_inputDO_%s_l2pen_%s_lr_%s'%(
        '-'.join(np.array(param['n_hidden']).astype(str)),
        param['n_input'],
        param['activation'].__name__,
        param['hidden_dropout'],
        param['input_dropout'],
        param['l2_penalty'],
        param['learning_rate'],
    )

    if figure_subfolder not in os.listdir(output_folder):
        os.mkdir(output_folder + figure_subfolder)
    return output_folder + figure_subfolder

def make_drug_folder(GDSC_drug_name, target_drug_name, figure_folder, target='TCGA'):
    folder = 'GDSC_%s_%s_%s'%(GDSC_drug_name, target, target_drug_name)
    if folder not in os.listdir(figure_folder):
        os.mkdir('%s/%s'%(figure_folder, folder))
    return folder


def make_network(param):
    # Makes succession of layers
    hidden_output = [[torch.nn.Dropout(param['hidden_dropout'] if i > 0 else param['input_dropout']),
                    torch.nn.Linear(e1, e2),
                    param['activation']()]
                 for i, (e1,e2) in enumerate(zip([param['n_input']] + param['n_hidden'][:-1], param['n_hidden']))]
    hidden_output = functools.reduce(lambda x,y: x+y, hidden_output)

    net = torch.nn.Sequential(
                *hidden_output,
                torch.nn.Dropout(param['hidden_dropout']),
                torch.nn.Linear(param['n_hidden'][-1], param['n_output'])
    )

    return net

def train_batch(x_batch, y_batch, device, net, loss_func, optimizer):
    
    x_batch = x_batch.to(device)
    y_batch = y_batch.to(device)
    
    net.train()
    optimizer.zero_grad()
    y_predict_train = net(x_batch)
    loss = loss_func(y_predict_train, y_batch)
    loss.backward()
    optimizer.step() 

    return loss
