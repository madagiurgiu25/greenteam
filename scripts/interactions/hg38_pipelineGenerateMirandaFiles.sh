#!/bin/bash

# bash /Users/Diana/Desktop/hubs/hg38_pipelineGenerateMirandaFiles.sh "/Users/Diana/Desktop/hubs" "/Users/Diana/Desktop/hubs/hg38" "/Users/Diana/Desktop/orthologs"


printf "\n\nCREATING INTERACTION AUXILLIARY FILES (hg38)\n\n"



SCRIPT_DIR_NAME=$1
FILE_DIR_NAME=$2
ORTHOLOGS_DIR=$3

printf "Your .py scripts are located in: $SCRIPT_DIR_NAME \n"
printf "Your input and output files are located in: $FILE_DIR_NAME\n\n"
printf "Your ortholog file is located in: $ORTHOLOGS_DIR\n\n"


printf "### Merge all miranda.json into 1 file: hg38_interactions_allDBs_and_pc.json\n"
cat $FILE_DIR_NAME/Mirbase_human_gencode_lncRNA_filtered_new.json > $FILE_DIR_NAME/hg38_interactions_allDBs_and_pc.json
cat $FILE_DIR_NAME/Mirbase_human_gencode_pc_filtered_new.json >> $FILE_DIR_NAME/hg38_interactions_allDBs_and_pc.json
cat $FILE_DIR_NAME/Mirbase_human_lncipedia_filtered.json >> $FILE_DIR_NAME/hg38_interactions_allDBs_and_pc.json
cat $FILE_DIR_NAME/Mirbase_human_noncode_filtered.json >> $FILE_DIR_NAME/hg38_interactions_allDBs_and_pc.json
sed -i -e 's/\]\[/,/g' $FILE_DIR_NAME/hg38_interactions_allDBs_and_pc.json  

printf "### Generate dictionary of genes / transcripts per each gene (unique ids) from all DBs: hg38_genesTranscripts_allDBs_and_pc.json\n"
python3 $SCRIPT_DIR_NAME/getListGenesTranscripts.py $FILE_DIR_NAME hg38_noncode_newkeys.json hg38_gencode_long_noncoding_newkeys.json hg38_lncipedia_newkeys.json hg38_protein_coding.json hg38_genesTranscripts_allDBs_and_pc.json

printf "### Create files of unique gene - transcript - mirna interactions (one interaction per line, incl. alignment information)\n" 
printf "... parsing gencode lncrna\n"
python3 $SCRIPT_DIR_NAME/parseMirandaRelDB.py $FILE_DIR_NAME/ Mirbase_human_gencode_lncRNA_filtered_new.json hg38_interactions_gencode_lncRNA.txt
printf "... parsing gencode pc\n"
python3 $SCRIPT_DIR_NAME/parseMirandaRelDB.py $FILE_DIR_NAME/ Mirbase_human_gencode_pc_filtered_new.json hg38_interactions_gencode_pc.txt
printf "... parsing lncipedia\n"
python3 $SCRIPT_DIR_NAME/parseMirandaRelDB.py $FILE_DIR_NAME/ Mirbase_human_lncipedia_filtered.json hg38_interactions_lncipedia.txt
printf "... parsing noncode\n"
python3 $SCRIPT_DIR_NAME/parseMirandaRelDB.py $FILE_DIR_NAME/ Mirbase_human_noncode_filtered.json hg38_interactions_noncode.txt

printf "### Merge all *_interactions_*.txt files into 1 file: hg38_interactionsAllGenes.txt. Will be read in MirandaRelDB module\n"
printf "Name_gene\tName_miRNA\tName_transcript\talign_score\tenergy\tmirna_start\tmirna_end\tlnc_start\tlnc_end\talign_len\tmirna_iden\tlncrna_iden\tmirna_alignment\talignment\tlncrna_alignment\n" > $FILE_DIR_NAME/hg38_interactionsAllGenes.txt
cat $FILE_DIR_NAME/hg38_interactions_gencode_lncRNA.txt >> $FILE_DIR_NAME/hg38_interactionsAllGenes.txt
cat $FILE_DIR_NAME/hg38_interactions_gencode_pc.txt >> $FILE_DIR_NAME/hg38_interactionsAllGenes.txt
cat $FILE_DIR_NAME/hg38_interactions_lncipedia.txt >> $FILE_DIR_NAME/hg38_interactionsAllGenes.txt
cat $FILE_DIR_NAME/hg38_interactions_noncode.txt >> $FILE_DIR_NAME/hg38_interactionsAllGenes.txt 

printf "### Count number of binding sites for each gene with each miRNA. Per line: gene - miRNA - number of binding sites\n"
printf "... counting gencode lncrna\n"
python3 $SCRIPT_DIR_NAME/getGeneMirnaInteractions.py $FILE_DIR_NAME/ Mirbase_human_gencode_lncRNA_filtered_new.json hg38_number_of_interactions_gencode_lncRNA.txt
printf "... counting gencode pc\n"
python3 $SCRIPT_DIR_NAME/getGeneMirnaInteractions.py $FILE_DIR_NAME/ Mirbase_human_gencode_pc_filtered_new.json hg38_number_of_interactions_gencode_pc.txt
printf "... counting lncipedia\n"
python3 $SCRIPT_DIR_NAME/getGeneMirnaInteractions.py $FILE_DIR_NAME/ Mirbase_human_lncipedia_filtered.json hg38_number_of_interactions_lncipedia.txt
printf "... counting noncode\n"
python3 $SCRIPT_DIR_NAME/getGeneMirnaInteractions.py $FILE_DIR_NAME/ Mirbase_human_noncode_filtered.json hg38_number_of_interactions_noncode.txt

printf "### Merge all *_number_of_interactions_*.txt files into 1 file: hg38_number_of_interactions.txt. Used for hub statistics\n"
printf "gene_id\tmiRNA\tbs_start\tbs_stop\tnr_interactions_for_bs_per_pair\tnr_unique_bs_per_pair\tnr_interactions_per_pair\ttotal_nr_unique_bs_per_gene\ttotal_nr_interactions_per_gene\n" > $FILE_DIR_NAME/hg38_number_of_interactions.txt
cat $FILE_DIR_NAME/hg38_number_of_interactions_gencode_lncRNA.txt >> $FILE_DIR_NAME/hg38_number_of_interactions.txt
cat $FILE_DIR_NAME/hg38_number_of_interactions_gencode_pc.txt >> $FILE_DIR_NAME/hg38_number_of_interactions.txt
cat $FILE_DIR_NAME/hg38_number_of_interactions_lncipedia.txt >> $FILE_DIR_NAME/hg38_number_of_interactions.txt
cat $FILE_DIR_NAME/hg38_number_of_interactions_noncode.txt >> $FILE_DIR_NAME/hg38_number_of_interactions.txt

printf "\n### Count number of gene targets for each miRNA: hg38_miRNA_number_of_targets.txt. Per line: miRNA - number of unique gene targets\n"
cat $FILE_DIR_NAME/hg38_number_of_interactions.txt| awk -F"\t" '{print $2}' | sort | uniq -c > $FILE_DIR_NAME/hg38_miRNA_number_of_targets.txt

printf "\n### Generate list of hub interactions: gene - miRNA - #interactions per pair\n"
printf "... non-strict hub interactions: hg38_hubInteractions_non_strict.txt\n"
printf "gene_id\tmiRNA\tnr_interactions_for_bs_per_pair\tnr_unique_bs_per_pair\tnr_interactions_per_pair\n" > $FILE_DIR_NAME/hg38_hubInteractions_non_strict.txt
tail -n+2 $FILE_DIR_NAME/hg38_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq >> $FILE_DIR_NAME/hg38_hubInteractions_non_strict.txt
printf "... non-strict hub interactions: hg38_hubInteractions_strict.txt\n"
printf "gene_id\tmiRNA\tnr_interactions_for_bs_per_pair\tnr_unique_bs_per_pair\tnr_interactions_per_pair\n" > $FILE_DIR_NAME/hg38_hubInteractions_strict.txt
tail -n+2 $FILE_DIR_NAME/hg38_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' && $5 > '1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq >> $FILE_DIR_NAME/hg38_hubInteractions_strict.txt

printf "\n### Generate list of hub genes with miRNA and their total number of binding sites and cummulated interactions\n"
printf "... non-strict hub list: hg38_hubGenes_total_interactions_bs_non_strict.txt\n"
tail -n+2 $FILE_DIR_NAME/hg38_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' {print $1"\t"$8"\t"$9}' | sort -n -k3 | uniq > $FILE_DIR_NAME/hg38_hubGenes_total_interactions_bs_non_strict.txt
printf "... strict hub list: hg38_hubGenes_total_interactions_bs_non_strict.txt\n"
tail -n+2 $FILE_DIR_NAME/hg38_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' && $5 > '1'  {print $1"\t"$8"\t"$9}' | sort -n -k3 | uniq > $FILE_DIR_NAME/hg38_hubGenes_total_interactions_bs_strict.txt

printf "\n### Generate list of gene features: hub_score (total number of interactions per hub), ortholog, list of mirna targets: hg38_geneFeatures.txt, hg38_geneFeatures.json\n"
python3 $SCRIPT_DIR/geneFeatures.py $FILE_DIR $ORTHOLOGS_DIR /Users/Diana/Desktop/overlap/ hg38_number_of_interactions.txt hg38_hubGenes_total_interactions_bs_non_strict.txt orthologs_genes_and_pc.txt hg38 hg38_geneFeatures.json

