import os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
import scipy
from copy import deepcopy
import uuid
from pickle import load, dump
import re
sns.set_style("whitegrid")
sns.set_context('paper')
from matplotlib import font_manager as fm, rcParams
fpath = os.path.join(rcParams["datapath"], "fonts/ttf/arial.ttf")
prop_label = fm.FontProperties(fname=fpath)
prop_label.set_size(30)
prop_ticks = fm.FontProperties(fname=fpath)
prop_ticks.set_size(25)
fname = os.path.split(fpath)[1]

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import Pipeline
from joblib import dump, load, Parallel, delayed
from statannot.statannot import add_stat_annotation
import torch
from skorch import NeuralNetClassifier, NeuralNetRegressor

sys.path.insert(0, '../read_data/')
from read_data import read_data
from read_GDSC_response import read_GDSC_response
from read_PDXE_response import read_PDXE_response
from reformat_df import reformat_df
import library_size_normalization

sys.path.insert(0, '../src/')
from clf_utils import make_network

from transact.pv_computation import PVComputation
from transact.interpolation import Interpolation
from transact.matrix_operations import _center_kernel, _right_center_kernel, _left_center_kernel
from transact.kernel_computer import KernelComputer
from transact.TRANSACT import TRANSACT