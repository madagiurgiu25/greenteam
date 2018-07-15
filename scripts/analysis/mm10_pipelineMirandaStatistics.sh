#!/bin/bash

printf "\n\nCREATING NECESSARY FILES\n\n"

# printf "\n### Merge all miranda.json into 1 file: mm10_interactions_allDBs_and_pc.json\n"
# cat Mirbase_mouse_gencode_lncRNA_filtered_new.json > mm10_interactions_allDBs_and_pc.json
# cat Mirbase_mouse_gencode_pc_filtered_new.json >> mm10_interactions_allDBs_and_pc.json
# cat Mirbase_mouse_noncode_filtered.json >> mm10_interactions_allDBs_and_pc.json
# sed -i -e 's/\]\[/,/g' mm10_interactions_allDBs_and_pc.json  

# printf "\n### Generate dictionary of genes / transcripts per each gene (unique ids) from all DBs: mm10_genesTranscripts_allDBs_and_pc.json\n"
# python3 ../getListGenesTranscripts.py ./ mm10_noncode_newkeys.json mm10_gencode_long_noncoding_newkeys.json mm10_protein_coding.json mm10_genesTranscripts_allDBs_and_pc.json

# printf "\n### Create files of unique gene - transcript - mirna interactions (one interaction per line, incl. alignment information)\n" 
# printf "... parsing gencode lncrna\n"
# python3 ../parseMirandaRelDB.py ./ Mirbase_mouse_gencode_lncRNA_filtered_new.json mm10_interactions_gencode_lncRNA.txt
# printf "... parsing gencode pc\n"
# python3 ../parseMirandaRelDB.py ./ Mirbase_mouse_gencode_pc_filtered_new.json mm10_interactions_gencode_pc.txt
# printf "... parsing noncode\n"
# python3 ../parseMirandaRelDB.py ./ Mirbase_mouse_noncode_filtered.json mm10_interactions_noncode.txt

# printf "\n### Merge all *_interactions_*.txt files into 1 file: mm10_interactionsAllGenes.txt. Will be read in MirandaRelDB module\n"
# printf "Name_gene\tName_miRNA\tName_transcript\talign_score\tenergy\tmirna_start\tmirna_end\tlnc_start\tlnc_end\talign_len\tmirna_iden\tlncrna_iden\tmirna_alignment\talignment\tlncrna_alignment\n" > mm10_interactionsAllGenes.txt
# cat mm10_interactions_gencode_lncRNA.txt >> mm10_interactionsAllGenes.txt
# cat mm10_interactions_gencode_pc.txt >> mm10_interactionsAllGenes.txt
# cat mm10_interactions_noncode.txt >> mm10_interactionsAllGenes.txt 

# printf "\n### Count number of binding sites for each gene with each miRNA. Per line: gene - miRNA - number of binding sites\n"
# printf "... counting gencode lncrna\n"
# python3 ../getGeneMirnaInteractions.py ./ Mirbase_mouse_gencode_lncRNA_filtered_new.json mm10_number_of_interactions_gencode_lncRNA.txt
# printf "... counting gencode pc\n"
# python3 ../getGeneMirnaInteractions.py ./ Mirbase_mouse_gencode_pc_filtered_new.json mm10_number_of_interactions_gencode_pc.txt
# printf "... counting noncode\n"
# python3 ../getGeneMirnaInteractions.py ./ Mirbase_mouse_noncode_filtered.json mm10_number_of_interactions_noncode.txt

# printf "\n### Merge all *_number_of_interactions_*.txt files into 1 file: mm10_number_of_interactions.txt. Used for hub statistics\n"
# cat mm10_number_of_interactions_gencode_lncRNA.txt > mm10_number_of_interactions.txt
# cat mm10_number_of_interactions_gencode_pc.txt >> mm10_number_of_interactions.txt
# cat mm10_number_of_interactions_noncode.txt >> mm10_number_of_interactions.txt

# printf "\n### Count number of gene targets for each miRNA: mm10_miRNA_number_of_targets.txt. Per line: miRNA - number of unique gene targets\n"
# cat mm10_number_of_interactions.txt| awk -F"\t" '{print $2}' | sort | uniq -c > mm10_miRNA_number_of_targets.txt

printf "\n\nCALCULATING STATISTICS\n\n"

printf "\n### Total number of interactions (filtered, but not processed)\n" 
printf "... total filtered predictions in gencode lncrna: "
cat Mirbase_mouse_gencode_lncRNA_filtered_new.tsv | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... total filtered predictions in gencode pc: "
cat Mirbase_mouse_gencode_pc_filtered_new.tsv | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... total filtered predictions in noncode: "
cat Mirbase_mouse_noncode_filtered.tsv | wc -l | awk -F"\t" '{print $1-'1'}'

printf "\n### Total number of interactions (filtered and processed - after mapping)\n"
printf "... total filtered predictions in gencode lncrna: "
cat mm10_interactions_gencode_lncRNA.txt | wc -l
printf "... total filtered predictions in gencode pc: "
cat mm10_interactions_gencode_pc.txt | wc -l
printf "... total filtered predictions in noncode: "
cat mm10_interactions_noncode.txt | wc -l

printf "\n### Total number of unique genes: "
cat mm10_number_of_interactions.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique LNC (noncode) genes: "
cat mm10_number_of_interactions_noncode.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique LNC (gencode) genes: "
cat mm10_number_of_interactions_gencode_lncRNA.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l
printf "... unique ENS (gencode) genes: "
cat mm10_number_of_interactions_gencode_pc.txt | awk -F"\t" '{print $1}' | sort | uniq | wc -l

printf "\n### Total number of unique miRNA: "
cat mm10_interactionsAllGenes.txt | awk -F"\t" '{print $2}' | sort | uniq | wc -l | awk -F"\t" '{print $1-'1'}'
printf "... unique miRNA with 1 gene target: "
cat mm10_miRNA_number_of_targets.txt | awk -F" " '$1 == '1' {print $2"\t"$1}' | sort -n -k2 | wc -l
printf "... unique miRNA with more than 1 gene target: "
cat mm10_miRNA_number_of_targets.txt | awk -F" " '$1 > '1' {print $2"\t"$1}' | sort -n -k2 | wc -l

printf "\n### Total miRNA-Gene interactions with 1 binding site (non-hub): "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 == '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | wc -l
printf "... unique LNC non-hub interactions: "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 == '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep LNC | wc -l
printf "... unique ENS non-hub interactions: "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 == '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep ENS | wc -l

printf "\n### Total miRNA-Gene interactions with more than 1 binding site at different genomic positions (hubs): "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | wc -l
printf "... unique LNC hub interactions: "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep LNC | wc -l
printf "... unique ENS hubs interactions: "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1"\t"$2"\t"$3 }' | sort -n -k3 | uniq | grep ENS | wc -l

printf "\n### Total number of unique hub genes: "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1 }' | sort -n -k3 | uniq | wc -l
printf "... unique LNC hub genes: "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1 }' | sort -n -k3 | uniq | grep LNC | wc -l
printf "... unique ENS hub genes: "
cat mm10_number_of_interactions.txt | sort -n -k3| awk -F"\t" '$3 > '1'{print $1 }' | sort -n -k3 | uniq | grep ENS | wc -l









# bash mm10_pipelineMirandaStatistics.sh




