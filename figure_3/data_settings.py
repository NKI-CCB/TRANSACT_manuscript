# Data sources
tissues = {
    'TCGA': ['All'],
    'GDSC': ['All']
}
projects = {
    'TCGA':[None],
    'GDSC': None
}

data_sources = ['GDSC', 'TCGA']

data_types = ['rnaseq']
genes_filtering = 'mini'
data_normalization = 'library_size'

source = 'GDSC'
target = 'TCGA'

# TRANSACT analysis
kernel_surname = 'rbf_gamma_0_0005'
kernel_name = 'rbf'
kernel_param = {
    'gamma': 0.0005
}

number_pc = {
    'source': 70,
    'target': 150
}
n_pv = 30
n_interpolation = 100
n_jobs = 20