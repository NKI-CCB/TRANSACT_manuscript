# Predicting clinical drug response from model systems by non-linear subspace-based transfer learning
## [Mourragui et al 2020]

Biorxiv pre-print : https://www.biorxiv.org/content/10.1101/2020.06.29.177139v2
<br/>
This folder contains the code necessary to:
- Download data (GDSC, PDXE and TCGA).
- Reproduce the main figures of the mananuscript (Fig 1,2,3,4 and 5) and the table (Table 2)

## Setting up

You need to have conda 4.8 or higher. Based on that, you can settup your environment using the following command:
<br/><br/>
<em> 
	conda create -n transact_figures python=3.6 <br/>
	conda activate transact_figures <br/>
	sh setup.sh <br/>
</em>

## Downloading the data

We provide a script to download the three public datasets used (GDSC, TCGA and PDXE). First you need to set-up the download environment:

<em>
	conda create -n data_download python=3.6 <br/>
	conda activate data_download <br/>
	cd ./download_data/ <br/>
	sh setup_download.sh <br/>
	sh download_data.sh <br/>
</em>
<br/>
HMF data can not be publicly shared. If you wish to reproduce our results, you must get in touch with the Hartwig Medical Foundation (https://www.hartwigmedicalfoundation.nl/en/) to get approval on data.

## Reproduce the results

Each folder contains the codes and notebook used to produce the figures. Each folder contains a sub-folder called './figures/' where all figures will be plotted in high resolution.
<br/><br/>
Of particular interest is the figure 4 that can be reproduced from subfolder "figure_4", but also contains the code to reproduce Table 1.

## Contact

In case you do not manage to reproduce the figures, you found a bug, or you have any question/suggestion on the code, please contact me either by GitHub, or by mail s [dot] mourragui [at] nki [dot] nl.