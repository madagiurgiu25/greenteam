import json
import pprint
import sys
import os
from pprint import pprint

'''
Generate dictionary of genes / transcripts per each gene (unique ids) from all DBs
input: *newkeys*.json files && [mm10_hg38]_protein_coding.json
output: mm10_genesTranscripts_allDBs_and_pc.json
'''
def loadJson(my_file):
	df_file = os.path.join(cwd, my_file)
	with open(df_file) as df_input:    
	    df = json.load(df_input)
	return df

def getGenesTranscripts(df):
	for gene in df:
		gene_id = gene["gene_id"]
		genes_transcripts[gene_id] = []
		for transcripts in gene["transcripts_list"]:
			transcript_id = transcripts["transcript_id"]
			genes_transcripts[gene_id].append(transcript_id)


if __name__ == '__main__':
    cwd = sys.argv[1]
    fin1 = sys.argv[2]
    fin2 = sys.argv[3]
    fin3 = sys.argv[4]
    if len(sys.argv) > 7:
        fin4 = sys.argv[5]
        fout = sys.argv[6]

        df1 = loadJson(fin1)
        df2 = loadJson(fin2)
        df3 = loadJson(fin3)
        df4 = loadJson(fin4)

        genes_transcripts = {}
        getGenesTranscripts(df1)
        getGenesTranscripts(df2)
        getGenesTranscripts(df3)
        getGenesTranscripts(df4)

        f_out = os.path.join(cwd, fout)
        output = open(f_out, "w")
        output.write(json.dumps(genes_transcripts, indent=4)) 
    else:
        fout = sys.argv[5]

        df1 = loadJson(fin1)
        df2 = loadJson(fin2)
        df3 = loadJson(fin3)

        genes_transcripts = {}
        getGenesTranscripts(df1)
        getGenesTranscripts(df2)
        getGenesTranscripts(df3)

        f_out = os.path.join(cwd, fout)
        output = open(f_out, "w")
        output.write(json.dumps(genes_transcripts, indent=4)) 

#python3 getListGenesTranscripts.py mm10/ mm10_noncode_newkeys.json mm10_gencode_lncRNA_newkeys.json mm10_lncrnadb_newkeys.json mm10_protein_coding.json mm10_genesTranscripts_allDBs_and_pc.json
#python3 getListGenesTranscripts.py hg38/ hg38_noncode_newkeys.json hg38_gencode_long_noncoding_newkeys.json hg38_lncipedia_newkeys.json hg38_protein_coding.json hg38_genesTranscripts_allDBs_and_pc.json
