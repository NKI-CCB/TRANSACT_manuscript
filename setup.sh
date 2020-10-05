# Activate Jupyter notebook
pip install ipykernel
python -m ipykernel install --user --name transact_figures --display-name "Python (TRANSACT_figures)"
conda install tzlocal

# Install kPRECISE
pip install transact_dr

# For notebooks
pip install matplotlib
pip install seaborn
conda install xarray
pip install umap-learn
pip install statannot
conda install netcdf4
pip install xlrd


# Install edgeR
conda install -c conda-forge r=3.5
conda install -c r rpy2
python install_edgeR.py