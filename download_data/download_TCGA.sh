eval "$(conda shell.bash hook)"
conda activate data_download

raw_folder='../data/TCGA/raw'
manifest_folder='./TCGA_manifests/'

mkdir -p $raw_folder
./gdc-client/bin/gdc-client 'download' -m $manifest_folder'manifest.txt' -d $raw_folder

# Supplement can be changed and retrieved from here https://academic.oup.com/bioinformatics/article/32/19/2891/2196464#supplementary-data
wget -O $raw_folder'response.zip' "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/32/19/10.1093_bioinformatics_btw344/4/btw344_supplementary_data.zip?Expires=1607528148&Signature=rqaiCbQpJSF2YCxht9jv6TdaFHrajCifzKO-CrHsC~qlTBwkvweF-usl1-v8~ihJPoCk0OgQD7aNpwPJiNnxka6LTBs4-C87DbEOlo2RARZwLg52rsB-aAXolty8S92QSY7DC3k3HrksQYpCI-b1xWfIfHWtSK6SEwIWnjvC6Hcpz9U0y6kyvk1P0Tvgkpt4q4eWGM3W0yFw68LnNjxh7xLNf9JRRfRKj2wPokR4mMt8BNoy15C7UVvblQY0gZkNy23Rlwi2pd0s2tQBF7~PD4yo3JpHrokd29Q~yf9dzwOvsxDxEz9NrGn8IfA~eshtvzvySlWFvOE21YZclkHpDg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"
unzip $raw_folder'response.zip' -d $raw_folder

python ./TCGA_dependencies/download_TCGA.py
python ./TCGA_dependencies/process_TCGA.py

rm -r $raw_folder

