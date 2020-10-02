mkdir ../data/GDSC/raw

wget -P ../data/GDSC/raw https://cog.sanger.ac.uk/cmp/download/rnaseq_20191101.zip
unzip ../data/GDSC/raw/rnaseq_20191101.zip -d ../data/GDSC/raw/
rm ../data/GDSC/raw/rnaseq_20191101.zip

wget -P ../data/GDSC https://cog.sanger.ac.uk/cmp/download/model_list_20191104.csv

eval "$(conda shell.bash hook)"
conda activate data_download
python ./GDSC_dependencies/process_GDSC.py

wget -P ../data/GDSC/response ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_25Feb20.xlsx
wget -P ../data/GDSC/response ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_25Feb20.xlsx

rm -r ../data/GDSC/raw