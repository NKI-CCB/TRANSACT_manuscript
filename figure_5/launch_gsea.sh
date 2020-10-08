## 2020/05/28 LAUNCH GSEA

IFS=","
input_file=$1
pathway="$2"
permutations=$3
while read -r drug_name drug_id
do
	echo $drug_name;
	echo $drug_id;
	echo '';
    sh gsea_one_gs.sh $pathway $drug_name $drug_id $permutations
done < $input_file