import json
import pprint
import sys
import os
from pprint import pprint
from collections import defaultdict

def getGenesMirna(df):
	for line in df:
		line_array = line.split("\t")
		gene_id = str(line_array[0])
		mirna_id = str(line_array[1])
		if gene_id in miranda_genes_mirna:
			miranda_genes_mirna[gene_id].append(mirna_id)
		else:
			miranda_genes_mirna[gene_id] = [mirna_id]

if __name__ == '__main__':
	fin_mirandaGeneList = sys.argv[1]
	fin_diffexpGeneList = sys.argv[2]
	fout = sys.argv[3]

	miranda_input = open(fin_mirandaGeneList, "r")
	diffexp_input = open(fin_diffexpGeneList, "r")
	output = open(fout, "w")

	miranda_genes_mirna = {}
	getGenesMirna(miranda_input)

	for line in diffexp_input:
		line_array = line.split("\t")
		gene = line_array[2]
		if gene in miranda_genes_mirna:
			mirnas = '|'.join(miranda_genes_mirna[gene])
			line = line.rstrip() + "\t" + 'yes' + "\t" + mirnas + "\n"
		else:
			line = line.rstrip() + "\t" + 'no' + "\n"
		output.write(line)

# python3 overlapMirandaDiffExp.py ./mm10/mm10_interactions_m103_all.txt ./mm10/genexp_signif_m103.tab ./mm10/genexp_signif_m103_overlap_miranda.tab
# python3 overlapMirandaDiffExp.py ./mm10/mm10_interactions_mlet7_all.txt ./mm10/genexp_signif_mlet7.tab ./mm10/genexp_signif_mlet7_overlap_miranda.tab
# python3 overlapMirandaDiffExp.py ./mm10/mm10_interactions_mlet7_unfiltered.txt ./mm10/genexp_signif_mlet7.tab ./mm10/genexp_signif_mlet7_overlap_miranda_unfiltered.tab
# python3 overlapMirandaDiffExp.py ./mm10/mm10_interactions_m103_unfiltered.txt ./mm10/genexp_signif_m103.tab ./mm10/genexp_signif_m103_overlap_miranda_unfiltered.tab
# python3 overlapMirandaDiffExp.py ./mm10/mm10_number_of_interactions.txt ./mm10/genexp_signif_m103.tab ./mm10/genexp_signif_m103_overlap_miranda_test.tab

