eval "$(conda shell.bash hook)"
conda activate data_download

raw_folder='../data/TCGA/raw'
manifest_folder='./TCGA_manifests/'

mkdir -p $raw_folder
./gdc-client/bin/gdc-client 'download' -m $manifest_folder'manifest.txt' -d $raw_folder

# Supplement can be changed and retrieved from here https://academic.oup.com/bioinformatics/article/32/19/2891/2196464#supplementary-data
wget -O $raw_folder'response.zip' "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/32/19/10.1093_bioinformatics_btw344/4/btw344_supplementary_data.zip?Expires=1597162292&Signature=C7PIMtzXJtaP3W5ke2Et0LfPq4IEweLqXeJlylOXY0E24BuHzDc2hGyIVx5L2JqnDnPz2LFuLzTNiWRRuRbh0OgZVGFH5SG3jM52vmVtQFNi59xKuZs1IqnsbGFRUz7L7ggoPxpK2uNKDCbmSRSzPOgIcGdTWZQyYHgCknx0E4X~FkEm4ucOWfL5dxbKIcloDsPn85eW8sul3r43x0g9Zt~ioJFG9ysoMen9jyxhFZiwlsq5GkupvKLKLcW1Mk6qZ2jXfBmCBZfP0m9ejYR6565na0AXpNK-iYx4xioLpfMWmK0Kda~lD7HLc0vZbVvCPgn9biok2L9CPisf8XHOhQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"
unzip $raw_folder'response.zip' -d $raw_folder

python ./TCGA_dependencies/download_TCGA.py
python ./TCGA_dependencies/process_TCGA.py

rm -r $raw_folder

