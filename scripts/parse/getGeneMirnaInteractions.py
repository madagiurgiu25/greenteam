import json
import pprint
import sys
import os
from pprint import pprint
from collections import defaultdict

def loadJson(my_file):
	gt_file = os.path.join(my_dir, my_file)
	with open(gt_file) as gt_input:    
	    df = json.load(gt_input)
	return df

def getGenesTranscripts(df):
    for gene in df:
        gene_id = gene['gene_id']
        for transcript in gene['transcript_list']:
            transcript_id = transcript['transcript_id']
            for mirna in transcript['interaction_list']:
                mirna_id = mirna['mirna']
                output.write(gene_id + "\t" + transcript_id + "\t" + mirna_id)

if __name__ == '__main__':
    fin = sys.argv[1]
    fout = sys.argv[2]

    my_dir = "./"
    df = loadJson(fin)

    output = open(fout, "w")
    getGenesTranscripts(df)

#python3 getGeneMirnaInteractions.py miranda_noncode_test.json noncode_gene_mirna_interaction_list.txt