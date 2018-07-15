import json
import pprint
import sys
import os
from pprint import pprint
from collections import defaultdict

def loadJson(my_file):
    df_file = os.path.join(cwd, my_file)
    with open(df_file) as df_input:    
        df = json.load(df_input)
    return df

def getGenesTranscripts(df):
    container = {}
    for gene in df:
        gene_id = gene['gene_id']
        found = []
        for transcript in gene['transcript_list']:
            transcript_id = transcript['transcript_id']
            for mirna in transcript['interaction_list']:
                mirna_id = mirna['mirna']
                for alignment in mirna['alignment_list']:
                    start = alignment['lnc_start']
                    end = alignment['lnc_end']
                    interaction = gene_id+"@@@"+mirna_id+"@@@"+start+"@@@"+end
                    if interaction in found:
                        #print("interaction found: " + interaction)
                        container[interaction] = container[interaction] + 1
                    else:
                        #print("new interaction: " + interaction)
                        container[interaction] = 1
                        found.append(interaction)

    for key, value in container.items():
        key_array = key.split("@@@")
        output.write(key_array[0] + "\t" + key_array[1] + "\t" + str(value) + "\n")

if __name__ == '__main__':
    cwd = sys.argv[1]
    fin = sys.argv[2]
    fout = sys.argv[3]

    df = loadJson(fin)

    f_out = os.path.join(cwd, fout)
    output = open(f_out, "w")
    getGenesTranscripts(df)

# python3 getGeneMirnaInteractions.py mm10/ mm10_interactions_allDBs_and_pc.json mm10_number_of_interactions_unique_genomic_positions.txt
# python3 getGeneMirnaInteractions.py mm10/ miranda_let7_unfiltered.json mm10_interactions_mlet7_unfiltered.txt
# python3 getGeneMirnaInteractions.py mm10/ miranda_miR103_unfiltered.json mm10_interactions_m103_unfiltered.txt

# python3 getGeneMirnaInteractions.py hg38/ Mirbase_human_gencode_lncRNA_filtered_new.json hg38_number_of_interactions_gencode_lncrna_test.txt



