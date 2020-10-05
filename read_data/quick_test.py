from read_data import read_data

# da: domain adaptation
tissues = {
    'PDXE': ['All'],
    'GDSC': ['All']
}
projects = {
    'PDXE':[None],
    'GDSC': None
}

data_sources = ['GDSC', 'PDXE']

data_types = ['fpkm']
genes_filtering = 'mini'
data_normalization = 'library_size' # Can be TPM, "library_size" or "log". Else will not have any influence.

source = 'GDSC'
target = 'PDXE'

folder = '../data/'
data_df = read_data(tissues=tissues,
                    data_types=[e for e in data_types],
                    projects=projects,
                    data_sources=data_sources,
                   folder_basis=folder)