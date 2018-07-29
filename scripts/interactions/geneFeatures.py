import json
import pprint
import sys
import os
import re
from pprint import pprint
from collections import defaultdict

gene_features = defaultdict(dict)

def getGenesMirna(df_file):
	df = open(df_file, "r")
	next(df)
	for line in df:
		line_array = line.replace("\n","").split("\t")
		gene_id = str(line_array[0]).split(".")[0] # ENSMUSG00000000197.8 --> ENSMUSG00000000197
		mirna_id = str(line_array[1])
		#gene_features[gene_id]['hub'] = None
		#gene_features[gene_id]['ortholog'] = None

		if gene_id in gene_features:
			if mirna_id not in gene_features[gene_id]['mirnas']:
				gene_features[gene_id]['mirnas'].append(mirna_id)
		else:
			gene_features[gene_id]['mirnas'] = []
			gene_features[gene_id]['mirnas'].append(mirna_id)

	for gene in gene_features.keys():
		gene_features[gene]['hub'] = None
		gene_features[gene]['ortholog'] = None
	df.close()

def getHubs(df_file):
	df = open(df_file, "r")
	for line in df:
		line_array = line.split("\t")
		gene_id = str(line_array[0]).split(".")[0]
		total_nr_interactions_per_gene = line_array[1]
		if gene_id in gene_features:
			gene_features[gene_id]['hub'] = int(total_nr_interactions_per_gene)
	df.close()

def getOrthologs(df_file, organism):
	df = open(df_file, "r")
	for line in df:
		line_array = line.split("\t")
		if organism == 'mm10':
			ortho_gene = line_array[0].rstrip()
			gene_id = line_array[1].rstrip()
		else:
			gene_id = line_array[0].rstrip()
			ortho_gene = line_array[1].rstrip()
		if gene_id in gene_features:
			gene_features[gene_id]['ortholog'] = ortho_gene
	df.close()

if __name__ == '__main__':
	hubs_dir 			= sys.argv[1]
	orthologs_dir		= sys.argv[2]
	overlap_dir			= sys.argv[3]
	fin_mirandaGeneList = sys.argv[4]
	fin_hubList 		= sys.argv[5]
	fin_orthoList		= sys.argv[6]
	organism			= sys.argv[7]
	fout 				= sys.argv[8]

	miranda_in 	= os.path.join(hubs_dir, fin_mirandaGeneList)
	hubs_in 	= os.path.join(hubs_dir, fin_hubList)
	ortho_in 	= os.path.join(orthologs_dir, fin_orthoList)
	out_json 	= os.path.join(overlap_dir, fout)
	out_txt		= os.path.join(overlap_dir, fout.replace(".json", ".txt"))
	
	output_json = open(out_json, "w")
	output_txt = open(out_txt, "w")
	output_txt.write("gene_id\tortholog\thub_score\tmirnas\n")

	getGenesMirna(miranda_in)
	getHubs(hubs_in)
	getOrthologs(ortho_in, organism)

	output_json.write(json.dumps(gene_features, indent=4)) 
	for gene in gene_features.keys():
		mirnas = '|'.join(gene_features[gene]['mirnas'])
		output_txt.write(gene+"\t"+str(gene_features[gene]['ortholog'])+"\t"+str(gene_features[gene]['hub'])+"\t"+mirnas+"\n")

# python3 geneFeatures.py /Users/Diana/Desktop/hubs/mm10/ /Users/Diana/Desktop/orthologs/ /Users/Diana/Desktop/overlap/ mm10_number_of_interactions.txt mm10_hubGenes_total_interactions_bs_non_strict.txt orthologs_genes_and_pc.txt mm10 mm10_geneFeatures.json
# python3 geneFeatures.py /Users/Diana/Desktop/hubs/hg38/ /Users/Diana/Desktop/orthologs/ /Users/Diana/Desktop/overlap/ hg38_number_of_interactions.txt hg38_hubGenes_total_interactions_bs_non_strict.txt orthologs_genes_and_pc.txt hg38 hg38_geneFeatures.json




