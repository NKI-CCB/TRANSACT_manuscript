import os, sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy
import re
from datetime import date
sns.set_style("whitegrid")
sns.set_context('paper')

from sklearn.metrics import roc_curve, plot_roc_curve, roc_auc_score
from statannot.statannot import add_stat_annotation

sys.path.insert(0, '../read_data/')
from read_data import read_data
from read_GDSC_response import read_GDSC_response
from read_TCGA_response import read_TCGA_response
from reformat_df import reformat_df
import library_size_normalization

sys.path.insert(0, '/home/s.mourragui/science/multi_omics_precise/src/TRANSACT/')
from transact.TRANSACT import TRANSACT

import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr

rpy2.robjects.numpy2ri.activate()
importr('pROC')