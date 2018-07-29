#!/bin/bash

# bash /Users/Diana/Desktop/mm10_pipelineGenerateMirandaStatistics.sh "/Users/Diana/Desktop/mm10"



printf "\n\nCALCULATING INTERACTION STATISTICS (mm10)\n\n"



FILE_DIR_NAME=$1
printf "Your input and output files are located in: $FILE_DIR_NAME\n\n"



printf "\n### Total number of interactions (filtered, but not processed)\n" 
printf "... total filtered predictions in gencode lncrna: "
tail -n+2  $FILE_DIR_NAME/Mirbase_mouse_gencode_lncRNA_filtered_new.tsv | wc -l | awk -F"\t" '{print $1}'
printf "... total filtered predictions in gencode pc: "
tail -n+2  $FILE_DIR_NAME/Mirbase_mouse_gencode_pc_filtered_new.tsv | wc -l | awk -F"\t" '{print $1}'
printf "... total filtered predictions in noncode: "
tail -n+2  $FILE_DIR_NAME/Mirbase_mouse_noncode_filtered.tsv | wc -l | awk -F"\t" '{print $1}'

printf "\n### Total number of interactions (filtered and processed - after mapping)\n"
printf "... total filtered predictions in gencode lncrna: "
cat $FILE_DIR_NAME/mm10_interactions_gencode_lncRNA.txt | wc -l
printf "... total filtered predictions in gencode pc: "
cat $FILE_DIR_NAME/mm10_interactions_gencode_pc.txt | wc -l
printf "... total filtered predictions in noncode: "
cat $FILE_DIR_NAME/mm10_interactions_noncode.txt | wc -l

printf "\n### Total number of unique genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique LNC (noncode) genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions_noncode.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique LNC (gencode) genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions_gencode_lncRNA.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique ENS (gencode) genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions_gencode_pc.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l

printf "\n### Total number of unique miRNA: "
cat $FILE_DIR_NAME/mm10_interactionsAllGenes.txt | awk -F"\t" '{print $2}' | sort | uniq | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... unique miRNA with 1 gene target: "
cat $FILE_DIR_NAME/mm10_miRNA_number_of_targets.txt | awk -F" " '$1 == '1' {print $2"\t"$1}' | sort -n -k2 | wc -l
printf "... unique miRNA with more than 1 gene target: "
cat $FILE_DIR_NAME/mm10_miRNA_number_of_targets.txt | awk -F" " '$1 > '1' {print $2"\t"$1}' | sort -n -k2 | wc -l

printf "\n### Total miRNA-gene interactions with 1 binding site (non-hub): "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 == '1'{print $1"\t"$2"\t"$6 }' | sort -n -k3 | uniq | wc -l
printf "... unique LNC non-hub interactions: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 == '1'{print $1"\t"$2"\t"$6 }' | sort -n -k3 | uniq | grep LNC | wc -l
printf "... unique ENS non-hub interactions: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 == '1'{print $1"\t"$2"\t"$6 }' | sort -n -k3 | uniq | uniq | grep ENS | wc -l

printf "\n### Total miRNA-gene interactions with number of unique binding sites per interaction pair > 1 (non-strict hubs): "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq | wc -l
printf "... unique LNC non-strict hub interactions: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq | grep LNC | wc -l
printf "... unique ENS non-strict hubs interactions: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq | grep ENS | wc -l 

printf "\n### Total miRNA-gene interactions with number of unique binding sites per interaction pair > 1 && number of interactions per binding site > 1 (strict hubs): "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' && $5 > '1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq | wc -l
printf "... unique LNC strict hub interactions: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' && $5 > '1'  {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq | grep LNC | wc -l
printf "... unique ENS strict hubs interactions: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' && $5 > '1'  {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq | grep ENS | wc -l 

printf "\n### Total number of unique non-hub genes: "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 == '1' {print $1}' | uniq | wc -l
printf "... unique LNC non-hub genes: "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 == '1' {print $1}' | uniq | grep LNC | wc -l
printf "... unique ENS non-hub genes: "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 == '1' {print $1}' | uniq | grep ENS | wc -l

printf "\n### Total number of unique non-strict hub genes: "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 > '1' {print $1}' | uniq | wc -l
printf "... unique LNC non-strict hub genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 > '1' {print $1}' | uniq | grep LNC | wc -l
printf "... unique ENS non-strict hub genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 > '1' {print $1}' | uniq | grep ENS | wc -l

printf "\n### Total number of unique strict hub genes: "
tail -n+2 $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 > '1' && $5 > '1' {print $1}' | uniq | wc -l
printf "... unique LNC strict hub genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 > '1' && $5 > '1' {print $1}' | uniq | grep LNC | wc -l
printf "... unique ENS strict hub genes: "
cat $FILE_DIR_NAME/mm10_number_of_interactions.txt | sort -n -k6| awk -F"\t" '$6 > '1' && $5 > '1' {print $1}' | uniq | grep ENS | wc -l


