import json
import pprint
import sys
import os
from pprint import pprint
from collections import defaultdict

'''
Create files of unique gene - transcript - mirna interactions (one interaction per line, incl. alignment information). Will be used for website (display miranda interactions as evidence)
input: miranda .json interaction file
output: mm10_interactions_*.txt
'''
def loadJson(my_file):
    df_file = os.path.join(cwd, my_file)
    with open(df_file) as df_input:    
        df = json.load(df_input)
    return df

def getGenesTranscripts(df):
    for gene in df:
        gene_id = gene['gene_id']
        for transcript in gene['transcript_list']:
            transcript_id = transcript['transcript_id']
            for mirna in transcript['interaction_list']:
                mirna_id = mirna['mirna']
                for align in mirna['alignment_list']:
                    align_score         = align['align_score']
                    energy              = align['energy'] 
                    mirna_start         = align['mirna_start']
                    mirna_end           = align['mirna_end']
                    lnc_start           = align['lnc_start']
                    lnc_end             = align['lnc_end']
                    align_len           = align['align_len']
                    mirna_iden          = align['mirna_iden']
                    lncrna_iden         = align['lncrna_iden']
                    mirna_alignment     = align['mirna_alignment']
                    alignment           = align['alignment']
                    lncrna_alignment    = align['lncrna_alignment']

                    output.write(gene_id + "\t" + mirna_id + "\t" + transcript_id 
                        + "\t" + align_score + "\t" + energy + "\t" + mirna_start + "\t" + mirna_end
                        + "\t" + lnc_start + "\t" + lnc_end + "\t" + align_len + "\t" + mirna_iden + lncrna_iden 
                        + "\t" + mirna_alignment + "\t" + alignment + "\t" + lncrna_alignment
                        + "\n")

if __name__ == '__main__':
    cwd = sys.argv[1]
    fin = sys.argv[2]
    fout = sys.argv[3]

    df = loadJson(fin)

    file_out = os.path.join(cwd, fout)
    output = open(file_out, "w")
    #output.write("Name_gene\tName_miRNA\tName_transcript\talign_score\tenergy\tmirna_start\tmirna_end\tlnc_start\tlnc_end"+"\talign_len\tmirna_iden\tlncrna_iden\tmirna_alignment\talignment\tlncrna_alignment\n")
    getGenesTranscripts(df)

#python3 parseMirandaRelDB.py ./ mm10_interactions_allDBs_and_pc.json mm10_interactionsAllGenes.txt

#python3 parseMirandaRelDB.py ./ Mirbase_mouse_gencode_lncRNA_filtered_new.json mm10_interactions_gencode_lncRNA.txt
#python3 parseMirandaRelDB.py ./ Mirbase_mouse_gencode_pc_filtered_new.json mm10_interactions_gencode_pc.txt
#python3 parseMirandaRelDB.py ./ Mirbase_mouse_noncode_filtered.json mm10_interactions_noncode.txt
