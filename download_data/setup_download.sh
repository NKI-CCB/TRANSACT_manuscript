# Set up the environment
conda create -y -n data_download python=3.6t
eval "$(conda shell.bash hook)"
conda activate data_download
conda install -y -c anaconda virtualenv
pip install xlrd==1.2.0
pip install -r ./TCGA_dependencies/requirements.txt

# Install GDC client
git clone https://github.com/NCI-GDC/gdc-client.git
cd gdc-client
pip install -r requirements.txt
python setup.py install
cd bin
./package
