import json
import pprint
import sys
import os
from pprint import pprint
from collections import defaultdict

my_dir = "./"

def loadJson(my_file):
	gt_file = os.path.join(my_dir, my_file)
	with open(gt_file) as gt_input:    
	    df = json.load(gt_input)
	return df

def getGenesTranscripts(df):
	for gene in df:
		gene_id = gene["gene_id"]
		genes_transcripts[gene_id] = []
		for transcripts in gene["transcripts_list"]:
			transcript_id = transcripts["transcript_id"]
			genes_transcripts[gene_id].append(transcript_id)
	#return genes_transcripts


if __name__ == '__main__':
    my_noncode_file = sys.argv[1]
    my_gencode_file = sys.argv[2]
    my_lncrnadb_file = sys.argv[3]

    df_noncode = loadJson(my_noncode_file)
    df_gencode = loadJson(my_gencode_file)
    df_lncrnadb = loadJson(my_lncrnadb_file)

    genes_transcripts = {}
    getGenesTranscripts(df_noncode)
    getGenesTranscripts(df_gencode)
    getGenesTranscripts(df_lncrnadb)

    #json_string = {}
    #json_string.append(gt_noncode)
    #json_string.append(gt_gencode)
    #json_string.append(gt_lncrnadb)

    out_file = os.path.join(my_dir, "genestranscripts.json")
    output = open(out_file, "w")
    output.write(json.dumps(genestranscripts, indent=4))
    #print(json.dumps(json_string, indent=4)) 

#python3 getListGenesTranscripts.py mm10_noncode_newkeys.json mm10_gencode_lncRNA_newkeys.json mm10_lncrnadb_newkeys.json