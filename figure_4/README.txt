FIGURE 4 & TABLE 1

This folder contains the scripts to reproduce:
- Figure 4
- Table 1
- Supp Fig 11
- Supp Fig 12


WARNINGS:
* If you have not yet set up the environment, please read "README.md" at the root of the
repository.
* If you have not downloaded the data, please report to the "download_data" subfolder.


In order to reproduce the aforementionned elements, please proceed as follows:

1) Launch the script launch_figure_4.sh. The scripts have been designed to store outputs
in ./figures/, but you may want to change that in the different scripts.
This script will already store in ./figures/TRANSACT/ the boxplots necessary for Fig 4B,
Fig 4C and Fig 4D.

2) For processing the different deep learning initializations (required for volcano
plot and table), load the Jupyter Notebook "instability_prediction_process". The only
parameter to change is the "data_type", which can be set either to "TCGA" or "HMF".
The different results will be stored in './figures/' and the cross-validated networks
will be taken out from '../deep_learning_CV/'. In case you have changed the folders 
beforehand, please do not forget to change folders in this notebook.

3) In order to process the p-values, load the Jupyter Notebook "process_p_values". As in 
"instability_prediction_process", you will need to change the parameter "data_type". This
notebook will output the Table 2 and the files required for Volcano plot.

4) In order to reproduce Fig 4A (volcano plot), load the Jupyter Notebook "volcano_plot".
This requires all three previous steps to be successful. Launching the notebook will
output the volcanoplot in ./figures.