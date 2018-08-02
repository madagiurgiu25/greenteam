#!/bin/bash

# bash /Users/Diana/Desktop/pipeline_overlaps.sh "/Users/Diana/Desktop" "/Users/Diana/Desktop/expression" "/Users/Diana/Desktop/interactions" "/Users/Diana/Desktop/orthologs" "/Users/Diana/Desktop/overlap"


printf "\n\n Overlapping diff exp - hubs - orthologs in human and mouse \n\n"


SCRIPT_DIR=$1
EXPRESSION_DIR=$2
HUBS_DIR=$3
ORTHOLOGS_DIR=$4
OVERLAP_DIR=$5

printf "Your .py scripts are located in: $SCRIPT_DIR \n"
printf "Your output files are located in: $OVERLAP_DIR \n\n"

printf "\n### For the gene expression, get non-strict hubs and miranda interaction overlaps (all diffexp, not only confidence set)\n"
printf "... mir-103 mouse all diffexp set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_mouse_103.txt mm10_geneFeatures.json mm10_mir103_all_overlapDiffExpInteractionsHubs.json
printf "... mir-103 human all diffexp set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_human_103.txt hg38_geneFeatures.json hg38_mir103_all_overlapDiffExpInteractionsHubs.json
printf "... mlet-7 mouse all diffexp set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "let-7" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_mouse_let7.txt mm10_geneFeatures.json mm10_let7_all_overlapDiffExpInteractionsHubs.json
printf "... hlet-7 human all diffexp set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "let-7" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_human_let7.txt hg38_geneFeatures.json hg38_let7_all_overlapDiffExpInteractionsHubs.json

printf "\n### For the gene expression, get non-strict hubs and miranda interaction overlaps (confident diffexp)\n"
printf "... mir-103 mouse confidence set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_mouse_103.txt_confident.txt mm10_geneFeatures.json mm10_mir103_confident_overlapDiffExpInteractionsHubs.json
printf "... mir-103 human confidence set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_human_103.txt_confident.txt hg38_geneFeatures.json hg38_mir103_confident_overlapDiffExpInteractionsHubs.json
printf "... mlet-7 mouse confidence set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "let-7" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_mouse_let7.txt_confident.txt mm10_geneFeatures.json mm10_let7_confident_overlapDiffExpInteractionsHubs.json
printf "... hlet-7 human confidence set\n"
python3 $SCRIPT_DIR/overlapGeneFeaturesDiffExp.py "let-7" $EXPRESSION_DIR $OVERLAP_DIR mirtrap_human_let7.txt_confident.txt hg38_geneFeatures.json hg38_let7_confident_overlapDiffExpInteractionsHubs.json

printf "\n### Number of differentially expressed genes with interactions"
printf "\n... mir-103 mouse confidence set:"
cat $OVERLAP_DIR/mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '$5=="TRUE" {print}' |wc -l
printf "\n... mir-103 human confidence set:"
cat $OVERLAP_DIR/hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '$5=="TRUE" {print}' |wc -l
printf "\n... mlet-7 mouse confidence set:"
cat $OVERLAP_DIR/mm10_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '$5=="TRUE" {print}' |wc -l
printf "\n... hlet-7 human confidence set:"
cat $OVERLAP_DIR/hg38_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '$5=="TRUE" {print}' |wc -l

printf "\n### Number of differentially expressed genes with interactions + are also hubs"
printf "\n... mir-103 mouse confidence set:"
cat $OVERLAP_DIR/mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |wc -l
cat $OVERLAP_DIR/mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print $1"\tmiR-103\t"$2"\t"$3"\t"$4"\t"$5}}'
printf "\n... mir-103 human confidence set:"
cat $OVERLAP_DIR/hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |wc -l
cat $OVERLAP_DIR/hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print $1"\tmiR-103\t"$2"\t"$3"\t"$4"\t"$5}}'
printf "\n... mlet-7 mouse confidence set:"
cat $OVERLAP_DIR/mm10_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |wc -l
cat $OVERLAP_DIR/mm10_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print $1"\tmiR-103\t"$2"\t"$3"\t"$4"\t"$5}}'
printf "\n... hlet-7 human confidence set:"
cat $OVERLAP_DIR/hg38_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |wc -l
cat $OVERLAP_DIR/hg38_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print $1"\tmiR-103\t"$2"\t"$3"\t"$4"\t"$5}}'

printf "\n### Number of differentially expressed genes with interactions + are also hubs + have interactions with mir-103 and let-7"
printf "\n... mir-103 mouse confidence set:"
cat $OVERLAP_DIR/mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |grep -i let-7 
printf "\n... mir-103 human confidence set:"
cat $OVERLAP_DIR/hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |grep -i let-7 
printf "\n... mlet-7 mouse confidence set:"
cat $OVERLAP_DIR/mm10_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |grep -i 'miR-103-\|miR-103[a-zA-Z]'
printf "\n... hlet-7 human confidence set:"
cat $OVERLAP_DIR/hg38_let7_confident_overlapDiffExpInteractionsHubs.txt | awk -F"\t" '{if($5=="TRUE" && $4!="None") {print}}' |grep -i 'miR-103-\|miR-103[a-zA-Z]'

# printf "\n### Put all orthologs into one file: orthologs_genes_and_pc.txt\n"
# cat $ORTHOLOGS_DIR/9606.10090.genelist.txt |awk -F"\t" '{print $1"\t"$3}' > $ORTHOLOGS_DIR/orthologs_genes_and_pc.txt
# cat $ORTHOLOGS_DIR/hg38_mm10_lncRNA_selected_independent.tsv |awk -F"\t" '{print $1"\t"$2}' >> $ORTHOLOGS_DIR/orthologs_genes_and_pc.txt
# sed -i -e 's/"//g' $ORTHOLOGS_DIR/orthologs_genes_and_pc.txt 

printf "\n### Calculate ortholog overlap\n"
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt mir103_h2m_confident2confident_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt mm10_mir103_all_overlapDiffExpInteractionsHubs.txt mir103_h2m_confident2all_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR hg38_mir103_all_overlapDiffExpInteractionsHubs.txt mm10_mir103_all_overlapDiffExpInteractionsHubs.txt mir103_h2m_all2all_orthoOverlap_allFeatures.txt

python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt mir103_m2h_confident2confident_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt hg38_mir103_all_overlapDiffExpInteractionsHubs.txt mir103_m2h_confident2all_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR mm10_mir103_all_overlapDiffExpInteractionsHubs.txt hg38_mir103_all_overlapDiffExpInteractionsHubs.txt mir103_m2h_all2all_orthoOverlap_allFeatures.txt

python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR hg38_let7_confident_overlapDiffExpInteractionsHubs.txt mm10_let7_confident_overlapDiffExpInteractionsHubs.txt let7_h2m_confident2confident_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR hg38_let7_confident_overlapDiffExpInteractionsHubs.txt mm10_let7_all_overlapDiffExpInteractionsHubs.txt let7_h2m_confident2all_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR hg38_let7_all_overlapDiffExpInteractionsHubs.txt mm10_let7_all_overlapDiffExpInteractionsHubs.txt let7_h2m_all2all_orthoOverlap_allFeatures.txt

python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR mm10_let7_confident_overlapDiffExpInteractionsHubs.txt hg38_let7_confident_overlapDiffExpInteractionsHubs.txt let7_m2h_confident2confident_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR mm10_let7_confident_overlapDiffExpInteractionsHubs.txt hg38_let7_all_overlapDiffExpInteractionsHubs.txt let7_m2h_confident2all_orthoOverlap_allFeatures.txt
python3 $SCRIPT_DIR/overlapOrtho.py $OVERLAP_DIR mm10_let7_all_overlapDiffExpInteractionsHubs.txt hg38_let7_all_overlapDiffExpInteractionsHubs.txt let7_m2h_all2all_orthoOverlap_allFeatures.txt

printf "\n### Diffexp orthologs + hubs\n"
printf "\n\n\n... mir103 h2m confident2confident:"
cat $OVERLAP_DIR/mir103_h2m_confident2confident_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/mir103_h2m_confident2confident_orthoOverlap_allFeatures.txt
printf "\n... mir103 h2m confident2confident + both hubs:"
cat $OVERLAP_DIR/mir103_h2m_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/mir103_h2m_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}'
printf "\n... mir103 h2m confident2all:"
cat $OVERLAP_DIR/mir103_h2m_confident2all_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/mir103_h2m_confident2all_orthoOverlap_allFeatures.txt
printf "\n... mir103 h2m confident2all + both hubs:"
cat $OVERLAP_DIR/mir103_h2m_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/mir103_h2m_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}'

printf "\n\n\n... mir103 m2h confident2confident:"
cat $OVERLAP_DIR/mir103_m2h_confident2confident_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/mir103_m2h_confident2confident_orthoOverlap_allFeatures.txt
printf "\n... mir103 m2h confident2confident + both hubs:"
cat $OVERLAP_DIR/mir103_m2h_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/mir103_m2h_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}'
printf "\n... mir103 m2h confident2all:"
cat $OVERLAP_DIR/mir103_m2h_confident2all_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/mir103_m2h_confident2all_orthoOverlap_allFeatures.txt
printf "\n... mir103 m2h confident2all + both hubs:"
cat $OVERLAP_DIR/mir103_m2h_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/mir103_m2h_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}'

printf "\n\n\n... let7 h2m confident2confident:"
cat $OVERLAP_DIR/let7_h2m_confident2confident_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/let7_h2m_confident2confident_orthoOverlap_allFeatures.txt 
printf "\n... let7 h2m confident2confident + both hubs:"
cat $OVERLAP_DIR/let7_h2m_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/let7_h2m_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}'
printf "\n... let7 h2m confident2all:"
cat $OVERLAP_DIR/let7_h2m_confident2all_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/let7_h2m_confident2all_orthoOverlap_allFeatures.txt 
printf "\n... let7 h2m confident2all + both hubs:"
cat $OVERLAP_DIR/let7_h2m_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/let7_h2m_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}' 

printf "\n\n\n... let7 m2h confident2confident:"
cat $OVERLAP_DIR/let7_m2h_confident2confident_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/let7_m2h_confident2confident_orthoOverlap_allFeatures.txt 
printf "\n... let7 m2h confident2confident + both hubs:"
cat $OVERLAP_DIR/let7_m2h_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/let7_m2h_confident2confident_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}'
printf "\n... let7 m2h confident2all:"
cat $OVERLAP_DIR/let7_m2h_confident2all_orthoOverlap_allFeatures.txt |wc -l
cat $OVERLAP_DIR/let7_m2h_confident2all_orthoOverlap_allFeatures.txt 
printf "\n... let7 m2h confident2all + both hubs:"
cat $OVERLAP_DIR/let7_m2h_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print}}' |wc -l
cat $OVERLAP_DIR/let7_m2h_confident2all_orthoOverlap_allFeatures.txt |awk -F"\t" '{if($4!="None" && $10!="None") {print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$10"\t"$11}}'











