import json
import pprint
import sys
import os
import re
from pprint import pprint
from collections import defaultdict

if __name__ == '__main__':
	target_mirna		= sys.argv[1]
	expression_dir 		= sys.argv[2]
	overlap_dir			= sys.argv[3]
	fin_diffexpGenes 	= sys.argv[4]
	fin_geneFeatures	= sys.argv[5]
	fout 				= sys.argv[6]

	diffexp_input 	= os.path.join(expression_dir, fin_diffexpGenes)
	genes_input 	= os.path.join(overlap_dir, fin_geneFeatures)
	out_json 		= os.path.join(overlap_dir, fout)
	out_txt			= os.path.join(overlap_dir, fout.replace(".json", ".txt"))

	gene_features 	= json.loads(open(genes_input,"r").read())
	diffexp_in 		= open(diffexp_input, "r")
	next(diffexp_in)
	output_json 	= open(out_json, "w")
	output_txt 		= open(out_txt, "w")
	output_txt.write("gene_id\tdiffexp_log2fc\tortholog\thub_score\tvalid_interaction\tmirnas\n")

	output = defaultdict(dict)
	for line in diffexp_in:
		line_array = line.split("\t")
		gene_id = line_array[0].split(".")[0]
		log2fc = line_array[2]
		output[gene_id]['log2fc'] = log2fc
		output[gene_id]['ortholog'] = None
		output[gene_id]['hub'] = 'None'
		output[gene_id]['valid_interaction'] = 'FALSE'
		output[gene_id]['mirnas'] = None

		#output_line = gene + "\t" + log2fc
		if gene_id in gene_features.keys():
			output[gene_id]['ortholog'] = gene_features[gene_id]['ortholog']
			if gene_features[gene_id]['hub'] == 'TRUE':
				output[gene_id]['hub'] = str(gene_features[gene_id]['total_nr_interactions_per_gene'])
			# check mirna
			interaction_check = 0
			for mirna in gene_features[gene_id]['mirnas']:
				if re.search(target_mirna, mirna):
					interaction_check = 1
					break
			if interaction_check == 1:
				output[gene_id]['valid_interaction'] = 'TRUE'
				output[gene_id]['mirnas'] = gene_features[gene_id]['mirnas']

	output_json.write(json.dumps(output, indent=4)) 
	output_json.close()
	for gene in output.keys():
		output_line = gene + "\t" + output[gene]['log2fc']
		if output[gene]['ortholog']:
			output_line = output_line +"\t"+output[gene]['ortholog']
		else:
			output_line = output_line +"\t" + 'None'
		output_line = output_line +"\t"+output[gene]['hub']+"\t"+output[gene]['valid_interaction']
		mirnas = None
		if output[gene]['mirnas']:
			mirnas = '|'.join(output[gene]['mirnas'])
			output_line = output_line +"\t" + mirnas
		else:
			output_line = output_line +"\t" + 'None'
		output_line = output_line +"\n"
		output_txt.write(output_line)

# python3 overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" expression/ overlap/ mirtrap_mouse_103.txt mm10_geneFeatures.json mm10_mir103_confident_overlapDiffExpInteractionsHubs.txt
# python3 overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" expression/ overlap/ mirtrap_human_103.txt hg38_geneFeatures.json hg38_mir103_confident_overlapDiffExpInteractionsHubs.txt
# python3 overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" expression/ overlap/ mirtrap_mouse_let7.txt mm10_geneFeatures.json mm10_let7_confident_overlapDiffExpInteractionsHubs.txt
# python3 overlapGeneFeaturesDiffExp.py "miR-103([a-zA-Z]?$|-)" expression/ overlap/ mirtrap_human_let7.txt hg38_geneFeatures.json hg38_let7_confident_overlapDiffExpInteractionsHubs.txt



