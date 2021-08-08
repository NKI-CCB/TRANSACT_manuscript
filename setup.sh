# Activate Jupyter notebook
conda create -n transact_figures python=3.6
conda activate transact_figures
pip install -r requirements.txt
pip install ipykernel
conda install rpy2
python -m ipykernel install --user --name transact_figures --display-name "Python (TRANSACT_figures)"
conda install tzlocal

# Install edgeR
conda install -c bioconda bioconductor-edger
# python install_edgeR.py

# Install TRANSACT
pip install transact_dr