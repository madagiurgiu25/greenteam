#!/bin/bash

# bash /Users/Diana/Desktop/hubs/mm10_pipelineGenerateMirandaFiles.sh "/Users/Diana/Desktop/hubs" "/Users/Diana/Desktop/hubs/mm10" "/Users/Diana/Desktop/orthologs"


printf "\n\nCREATING INTERACTION AUXILLIARY FILES (mm10)\n\n"



SCRIPT_DIR=$1
FILE_DIR=$2
ORTHOLOGS_DIR=$3

printf "Your .py scripts are located in: $SCRIPT_DIR \n"
printf "Your input and output files are located in: $FILE_DIR\n\n"
printf "Your ortholog file is located in: $ORTHOLOGS_DIR\n\n"




printf "\n### Merge all miranda.json into 1 file: mm10_interactions_allDBs_and_pc.json\n"
cat $FILE_DIR/Mirbase_mouse_gencode_lncRNA_filtered_new.json > $FILE_DIR/mm10_interactions_allDBs_and_pc.json
cat $FILE_DIR/Mirbase_mouse_gencode_pc_filtered_new.json >> $FILE_DIR/mm10_interactions_allDBs_and_pc.json
cat $FILE_DIR/Mirbase_mouse_noncode_filtered.json >> $FILE_DIR/mm10_interactions_allDBs_and_pc.json
sed -i -e 's/\]\[/,/g' $FILE_DIR/mm10_interactions_allDBs_and_pc.json  

printf "\n### Generate dictionary of genes / transcripts per each gene (unique ids) from all DBs: mm10_genesTranscripts_allDBs_and_pc.json\n"
python3 $SCRIPT_DIR/getListGenesTranscripts.py $FILE_DIR mm10_noncode_newkeys.json mm10_gencode_lncRNA_newkeys.json mm10_protein_coding.json mm10_genesTranscripts_allDBs_and_pc_test.json

printf "\n### Create files of unique gene - transcript - mirna interactions (one interaction per line, incl. alignment information)\n" 
printf "... parsing gencode lncrna\n"
python3 $SCRIPT_DIR/parseMirandaRelDB.py $FILE_DIR Mirbase_mouse_gencode_lncRNA_filtered_new.json mm10_interactions_gencode_lncRNA.txt
printf "... parsing gencode pc\n"
python3 $SCRIPT_DIR/parseMirandaRelDB.py $FILE_DIR Mirbase_mouse_gencode_pc_filtered_new.json mm10_interactions_gencode_pc.txt
printf "... parsing noncode\n"
python3 $SCRIPT_DIR/parseMirandaRelDB.py $FILE_DIR Mirbase_mouse_noncode_filtered.json mm10_interactions_noncode.txt

printf "\n### Merge all *_interactions_*.txt files into 1 file: mm10_interactionsAllGenes.txt. Will be read in MirandaRelDB module\n"
printf "Name_gene\tName_miRNA\tName_transcript\talign_score\tenergy\tmirna_start\tmirna_end\tlnc_start\tlnc_end\talign_len\tmirna_iden\tlncrna_iden\tmirna_alignment\talignment\tlncrna_alignment\n" > $FILE_DIR/mm10_interactionsAllGenes.txt
cat $FILE_DIR/mm10_interactions_gencode_lncRNA.txt >> $FILE_DIR/mm10_interactionsAllGenes.txt
cat $FILE_DIR/mm10_interactions_gencode_pc.txt >> $FILE_DIR/mm10_interactionsAllGenes.txt
cat $FILE_DIR/mm10_interactions_noncode.txt >> $FILE_DIR/mm10_interactionsAllGenes.txt 

printf "\n### Count number of binding sites for each gene with each miRNA. Per line: (1)-gene,(2)-miRNA,(3)-binding site start,(4)-binding site stop,(5)-number of interaction for current binding site for pair,(6)-number of unique binding sites per pair,(7)-number of interactions per pair,(8)-total number of unique binding sites for gene,(9)-total number of interactions per gene\n"
printf "... counting gencode lncrna\n"
python3 $SCRIPT_DIR/getGeneMirnaInteractions.py $FILE_DIR Mirbase_mouse_gencode_lncRNA_filtered_new.json mm10_number_of_interactions_gencode_lncRNA.txt
printf "... counting gencode pc\n"
python3 $SCRIPT_DIR/getGeneMirnaInteractions.py $FILE_DIR Mirbase_mouse_gencode_pc_filtered_new.json mm10_number_of_interactions_gencode_pc.txt
printf "... counting noncode\n"
python3 $SCRIPT_DIR/getGeneMirnaInteractions.py $FILE_DIR Mirbase_mouse_noncode_filtered.json mm10_number_of_interactions_noncode.txt

printf "\n### Merge all *_number_of_interactions_*.txt files into 1 file: mm10_number_of_interactions.txt. Used for hub statistics\n"
printf "gene_id\tmiRNA\tbs_start\tbs_stop\tnr_interactions_for_bs_per_pair\tnr_unique_bs_per_pair\tnr_interactions_per_pair\ttotal_nr_unique_bs_per_gene\ttotal_nr_interactions_per_gene\n" > $FILE_DIR/mm10_number_of_interactions.txt
cat $FILE_DIR/mm10_number_of_interactions_gencode_lncRNA.txt >> $FILE_DIR/mm10_number_of_interactions.txt
cat $FILE_DIR/mm10_number_of_interactions_gencode_pc.txt >> $FILE_DIR/mm10_number_of_interactions.txt
cat $FILE_DIR/mm10_number_of_interactions_noncode.txt >> $FILE_DIR/mm10_number_of_interactions.txt

printf "\n### Count number of gene targets for each miRNA: mm10_miRNA_number_of_targets.txt. Per line: miRNA - number of unique gene targets\n"
cat $FILE_DIR/mm10_number_of_interactions.txt| awk -F"\t" '{print $2}' | sort | uniq -c > $FILE_DIR/mm10_miRNA_number_of_targets.txt

printf "\n### Generate list of hub interactions: gene - miRNA - #interactions per pair\n"
printf "... non-strict hub interactions: mm10_hubInteractions_non_strict.txt\n"
printf "gene_id\tmiRNA\tnr_interactions_for_bs_per_pair\tnr_unique_bs_per_pair\tnr_interactions_per_pair\n" > $FILE_DIR/mm10_hubInteractions_non_strict.txt
tail -n+2 $FILE_DIR/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq >> $FILE_DIR/mm10_hubInteractions_non_strict.txt
printf "... non-strict hub interactions: mm10_hubInteractions_strict.txt\n"
printf "gene_id\tmiRNA\tnr_interactions_for_bs_per_pair\tnr_unique_bs_per_pair\tnr_interactions_per_pair\n" > $FILE_DIR/mm10_hubInteractions_strict.txt
tail -n+2 $FILE_DIR/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' && $5 > '1' {print $1"\t"$2"\t"$5"\t"$6"\t"$7}' | sort -n -k4 | uniq >> $FILE_DIR/mm10_hubInteractions_strict.txt

printf "\n### Generate list of hub genes with miRNA and their total number of binding sites and cummulated interactions\n"
printf "... non-strict hub list: mm10_hubGenes_total_interactions_bs_non_strict.txt\n"
tail -n+2 $FILE_DIR/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' {print $1"\t"$8"\t"$9}' | sort -n -k3 | uniq > $FILE_DIR/mm10_hubGenes_total_interactions_bs_non_strict.txt
printf "... strict hub list: mm10_hubGenes_total_interactions_bs_non_strict.txt\n"
tail -n+2 $FILE_DIR/mm10_number_of_interactions.txt | sort -n -k6 | awk -F"\t" '$6 >'1' && $5 > '1'  {print $1"\t"$8"\t"$9}' | sort -n -k3 | uniq > $FILE_DIR/mm10_hubGenes_total_interactions_bs_strict.txt

printf "\n### Generate list of gene features: hub_score (total number of interactions per hub), ortholog, list of mirna targets: mm10_geneFeatures.txt, mm10_geneFeatures.json\n"
python3 $SCRIPT_DIR/geneFeatures.py $FILE_DIR $ORTHOLOGS_DIR /Users/Diana/Desktop/overlap/ mm10_number_of_interactions.txt mm10_hubGenes_total_interactions_bs_non_strict.txt orthologs_genes_and_pc.txt mm10 mm10_geneFeatures.json

