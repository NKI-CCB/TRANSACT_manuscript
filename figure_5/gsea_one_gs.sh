## GSEA ONE GENE SET
# 
# Launch GSEA for one drug and one gene set.

# cd ./GSEA_4.0.3;
pathway="$1"
drug="$2"
drug_id="$3"
permutations=$4

output_folder="./gsea_output/"$pathway"_"$permutations"_perm"
mkdir $output_folder;
./GSEA_Linux_4.0.3/gsea-cli.sh GSEAPreranked \
						-gmx 'ftp.broadinstitute.org://pub/gsea/gene_sets/'$pathway'.v7.1.symbols.gmt' \
						-collapse No_Collapse\
						-mode Max_probe\
						-norm meandiv\
						-nperm $4\
						-rnk './gsea_input/GDSC_TCGA_'$drug'_'$drug_id'_loading.rnk'\
						-scoring_scheme weighted_p2\
						-rpt_label 'drug_'$drug'_id_'$drug_id'_pathway_'$pathway\
						-create_svgs false\
						-include_only_symbols true\
						-make_sets true\
						-plot_top_x 20\
						-rnd_seed timestamp\
						-set_max 500\
						-set_min 15\
						-zip_report false\
						-out $output_folder;