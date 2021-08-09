# Normalization
with_mean = True
with_std = True

# domain adaptation
tissues = {
    'PDXE': ['All'],
    'GDSC': ['All']
}
projects = {
    'PDXE':None,
    'GDSC': None
}

data_sources = ['GDSC', 'PDXE']

data_types = ['fpkm']
genes_filtering = 'mini'
data_normalization = 'library_size'

source = 'GDSC'
target = 'PDXE'

# Folder where CV has been saved
output_combat_cv_folder = ''
output_uncorrected_cv_folder = ''
random_state = 183627362