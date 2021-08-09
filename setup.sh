# Activate Jupyter notebook
conda create -n transact_figures python=3.6
conda activate transact_figures
pip install -r requirements.txt
pip install ipykernel
conda install -c r r=3.5.1
conda install -c r rpy2
python -m ipykernel install --user --name transact_figures --display-name "Python (TRANSACT_figures)"
conda install tzlocal
conda install -c conda-forge umap-learn
pip install seaborn scikit-learn statannot torch skorch

# Install edgeR
# conda install -c bioconda bioconductor-edger
python install_edgeR.py

# Install TRANSACT
pip install transact_dr