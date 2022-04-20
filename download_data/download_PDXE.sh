mkdir -p ../data/PDXE/raw

wget -P ../data/PDXE/raw https://static-content.springer.com/esm/art%3A10.1038%2Fnm.3954/MediaObjects/41591_2015_BFnm3954_MOESM10_ESM.xlsx

eval "$(conda shell.bash hook)"
conda activate data_download
python ./PDXE_dependencies/download_PDXE.py
python ./PDXE_dependencies/process_PDXE.py

rm -r ../data/PDXE/raw
