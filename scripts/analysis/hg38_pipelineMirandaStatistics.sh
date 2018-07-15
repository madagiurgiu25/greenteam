#!/bin/bash

printf "\n\nCREATING NECESSARY FILES\n\n"

# printf "### Merge all miranda.json into 1 file: hg38_interactions_allDBs_and_pc.json\n"
# cat Mirbase_human_gencode_lncRNA_filtered_new.json > hg38_interactions_allDBs_and_pc.json
# cat Mirbase_human_gencode_pc_filtered_new.json >> hg38_interactions_allDBs_and_pc.json
# cat Mirbase_human_lncipedia_filtered.json >> hg38_interactions_allDBs_and_pc.json
# cat Mirbase_human_noncode_filtered.json >> hg38_interactions_allDBs_and_pc.json
# sed -i -e 's/\]\[/,/g' hg38_interactions_allDBs_and_pc.json  

# printf "### Generate dictionary of genes / transcripts per each gene (unique ids) from all DBs: hg38_genesTranscripts_allDBs_and_pc.json\n"
# python3 ../getListGenesTranscripts.py ./ hg38_noncode_newkeys.json hg38_gencode_long_noncoding_newkeys.json hg38_lncipedia_newkeys.json hg38_protein_coding.json hg38_genesTranscripts_allDBs_and_pc.json

# printf "### Create files of unique gene - transcript - mirna interactions (one interaction per line, incl. alignment information)\n" 
# printf "... parsing gencode lncrna\n"
# python3 ../parseMirandaRelDB.py ./ Mirbase_human_gencode_lncRNA_filtered_new.json hg38_interactions_gencode_lncRNA.txt
# printf "... parsing gencode pc\n"
# python3 ../parseMirandaRelDB.py ./ Mirbase_human_gencode_pc_filtered_new.json hg38_interactions_gencode_pc.txt
# printf "... parsing lncipedia\n"
# python3 ../parseMirandaRelDB.py ./ Mirbase_human_lncipedia_filtered.json hg38_interactions_lncipedia.txt
# printf "... parsing noncode\n"
# python3 ../parseMirandaRelDB.py ./ Mirbase_human_noncode_filtered.json hg38_interactions_noncode.txt

# printf "### Merge all *_interactions_*.txt files into 1 file: hg38_interactionsAllGenes.txt. Will be read in MirandaRelDB module\n"
# printf "Name_gene\tName_miRNA\tName_transcript\talign_score\tenergy\tmirna_start\tmirna_end\tlnc_start\tlnc_end\talign_len\tmirna_iden\tlncrna_iden\tmirna_alignment\talignment\tlncrna_alignment\n" > hg38_interactionsAllGenes.txt
# cat hg38_interactions_gencode_lncRNA.txt >> hg38_interactionsAllGenes.txt
# cat hg38_interactions_gencode_pc.txt >> hg38_interactionsAllGenes.txt
# cat hg38_interactions_lncipedia.txt >> hg38_interactionsAllGenes.txt
# cat hg38_interactions_noncode.txt >> hg38_interactionsAllGenes.txt 

printf "### Count number of binding sites for each gene with each miRNA. Per line: gene - miRNA - number of binding sites\n"
printf "... counting gencode lncrna\n"
python3 ../getGeneMirnaInteractions.py ./ Mirbase_human_gencode_lncRNA_filtered_new.json hg38_number_of_interactions_gencode_lncRNA.txt
printf "... counting gencode pc\n"
python3 ../getGeneMirnaInteractions.py ./ Mirbase_human_gencode_pc_filtered_new.json hg38_number_of_interactions_gencode_pc.txt
printf "... counting lncipedia\n"
python3 ../getGeneMirnaInteractions.py ./ Mirbase_human_lncipedia_filtered.json hg38_number_of_interactions_lncipedia.txt
printf "... counting noncode\n"
python3 ../getGeneMirnaInteractions.py ./ Mirbase_human_noncode_filtered.json hg38_number_of_interactions_noncode.txt

printf "### Merge all *_number_of_interactions_*.txt files into 1 file: hg38_number_of_interactions.txt. Used for hub statistics\n"
cat hg38_number_of_interactions_gencode_lncRNA.txt > hg38_number_of_interactions.txt
cat hg38_number_of_interactions_gencode_pc.txt >> hg38_number_of_interactions.txt
cat hg38_number_of_interactions_lncipedia.txt >> hg38_number_of_interactions.txt
cat hg38_number_of_interactions_noncode.txt >> hg38_number_of_interactions.txt

printf "\n### Count number of gene targets for each miRNA: hg38_miRNA_number_of_targets.txt. Per line: miRNA - number of unique gene targets\n"
cat hg38_number_of_interactions.txt| awk -F"\t" '{print $2}' | sort | uniq -c > hg38_miRNA_number_of_targets.txt

printf "\n\nCALCULATING STATISTICS\n\n"

printf "\n### Total number of interactions (filtered, but not processed)\n" 
printf "... total filtered predictions in gencode lncrna: "
cat Mirbase_human_gencode_lncRNA_filtered_new.tsv | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... total filtered predictions in gencode pc: "
cat Mirbase_human_gencode_pc_filtered_new.tsv | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... total filtered predictions in lncipedia: "
cat Mirbase_human_lncipedia_filtered.tsv | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... total filtered predictions in noncode: "
cat Mirbase_human_noncode_filtered.tsv | wc -l | awk -F"\t" '{print $1-'1'}'

printf "\n### Total number of interactions (filtered and processed - after mapping)\n"
printf "... total filtered predictions in gencode lncrna: "
cat hg38_interactions_gencode_lncRNA.txt | wc -l
printf "... total filtered predictions in gencode pc: "
cat hg38_interactions_gencode_pc.txt | wc -l
printf "... total filtered predictions in lncipedia: "
cat hg38_interactions_lncipedia.txt | wc -l
printf "... total filtered predictions in noncode: "
cat hg38_interactions_noncode.txt | wc -l

printf "\n### Total number of unique genes: "
cat hg38_number_of_interactions.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique LNC (lncipedia) genes: "
cat hg38_number_of_interactions_lncipedia.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique LNC (noncode) genes: "
cat hg38_number_of_interactions_noncode.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique LNC (gencode) genes: "
cat hg38_number_of_interactions_gencode_lncRNA.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique ENS (gencode) genes: "
cat hg38_number_of_interactions_gencode_pc.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l

printf "\n### Total number of unique miRNA: "
cat hg38_interactionsAllGenes.txt | awk -F"\t" '{print $2}' | sort | uniq | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... unique miRNA with 1 gene target: "
cat hg38_miRNA_number_of_targets.txt | awk -F" " '$1 == '1' {print $2"\t"$1}' | sort -n -k2 | wc -l
printf "... unique miRNA with more than 1 gene target: "
cat hg38_miRNA_number_of_targets.txt | awk -F" " '$1 > '1' {print $2"\t"$1}' | sort -n -k2 | wc -l

printf "\n### Total miRNA-Gene interactions with 1 binding site (non-hub): "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 == '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | wc -l
printf "... unique LNC non-hub interactions: "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 == '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep LNC | wc -l
printf "... unique ENS non-hub interactions: "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 == '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep ENS | wc -l

printf "\n### Total miRNA-Gene interactions with more than 1 binding site at different genomic positions (hubs): "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | wc -l
printf "... unique LNC hub interactions: "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep LNC | wc -l
printf "... unique ENS hubs interactions: "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep ENS | wc -l

printf "\n### Total number of unique hub genes: "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1 }' | sort -n -k3 | uniq | wc -l
printf "... unique LNC hub genes: "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1 }' | sort -n -k3 | uniq | grep LNC | wc -l
printf "... unique ENS hub genes: "
cat hg38_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1 }' | sort -n -k3 | uniq | grep ENS | wc -l



# bash hg38_pipelineMirandaStatistics.sh




