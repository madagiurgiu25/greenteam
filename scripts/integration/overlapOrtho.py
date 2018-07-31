import json
import pprint
import sys
import os
import re
from pprint import pprint	

def calculateOrthologs(df_query_file, df_target_file):
	df_query = open(df_query_file, "r")
	next(df_query)

	for query in df_query:
		query_array = query.split("\t")
		query_gene = query_array[0]
		target_gene = query_array[2]

		output_line = query

		if target_gene != "None":
			df_target = open(df_target_file, "r")
			next(df_target)
			for target in df_target:
				target_array = target.split("\t")
				gene = target_array[0]
				if gene == target_gene:
					print("YUHUUUU: " + query_gene + " with " + target_gene + " and both are diff exp")
					output_line = output_line.rstrip().replace("\n", "") + "\t" + target.rstrip().replace("\n", "") + "\n"
					output.write(output_line)
			df_target.close()
	df_query.close()
	

if __name__ == '__main__':
	overlap_dir		= sys.argv[1]
	fin_query 		= sys.argv[2]
	fin_target 		= sys.argv[3]
	fout 			= sys.argv[4]

	query_in 		= os.path.join(overlap_dir, fin_query)
	target_in 		= os.path.join(overlap_dir, fin_target)
	out 			= os.path.join(overlap_dir, fout)

	output 			= open(out, "w")	

	calculateOrthologs(query_in, target_in)

# python3 overlapOrtho.py overlap/ hg38_let7_confident_overlapDiffExpInteractionsHubs.txt mm10_let7_confident_overlapDiffExpInteractionsHubs.txt let7_h2m_confident2confident_orthoOverlap_allFeatures.txt
# python3 overlapOrtho.py overlap/ hg38_let7_confident_overlapDiffExpInteractionsHubs.txt mm10_let7_all_overlapDiffExpInteractionsHubs.txt let7_h2m_confident2all_orthoOverlap_allFeatures.txt