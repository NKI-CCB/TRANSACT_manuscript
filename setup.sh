# Activate Jupyter notebook
pip install -r requirements.txt
pip install ipykernel
python -m ipykernel install --user --name transact_figures --display-name "Python (TRANSACT_figures)"
conda install tzlocal

# Install TRANSACT
pip install transact_dr