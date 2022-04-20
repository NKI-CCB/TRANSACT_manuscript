# Set up the environment
conda create -n data_download python=3.6
eval "$(conda shell.bash hook)"
conda activate data_download
conda install -c anaconda virtualenv
conda install setuptools
pip install xlrd==1.2.0
pip install -r ./TCGA_dependencies/requirements.txt

# Install GDC client
git clone https://github.com/NCI-GDC/gdc-client.git
cd gdc-client
pip install -r requirements.txt
python setup.py install
cd bin
./package
