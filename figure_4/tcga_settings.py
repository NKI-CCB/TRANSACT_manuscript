data_type = 'TCGA'
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

source = 'GDSC'
target = 'TCGA'

with_mean = True
with_std = True

kernel_surname = 'rbf_gamma_0_0005'
kernel_name = 'rbf'
kernel_param = {
    'gamma': 0.0005
}

ROC_combat_dl_ci = {}
ROC_transact_ci = {}
ROC_precise_ci = {}
ROC_dl_ci = {}

ROC_combat_dl_p_val = {}
ROC_transact_p_val = {}
ROC_precise_p_val = {}
ROC_dl_p_val = {}

response_size = {}

number_pc = {
    'source': 70,
    'target': 150
}

n_pv = 30
n_interpolation = 100

test = 'Mann-Whitney-ls'