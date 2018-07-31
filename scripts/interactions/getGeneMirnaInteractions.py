import json
import pprint
import sys
import os
from pprint import pprint
from collections import defaultdict

'''
Count number of binding sites for each gene with each miRNA. 
Per line: (1)-gene,(2)-miRNA,(3)-binding site start,(4)-binding site stop,(5)-number of interaction for current binding site for pair,
(6)-number of unique binding sites per pair,(7)-number of interactions per pair,
(8)-total number of unique binding sites for gene,(9)-total number of interactions per gene

input: miranda.json
output: mm10_number_of_interactions_*.txt
'''

def loadJson(my_file):
    df_file = os.path.join(cwd, my_file)
    with open(df_file) as df_input:    
        df = json.load(df_input)
    return df

def getInteractions(df):
    for gene in df:
        gene_id = gene['gene_id']
        nesteddict_gene_total[gene_id]['total_unique_bs'] = 0
        nesteddict_gene_total[gene_id]['total_nr_interactions'] = 0
        found_interaction_pair = []
        found_bs = []
        for transcript in gene['transcript_list']:
            transcript_id = transcript['transcript_id']
            for mirna in transcript['interaction_list']:
                mirna_id = mirna['mirna']
                interaction_pair = gene_id + "@@@" + mirna_id
                if interaction_pair in found_interaction_pair:
                    for alignment in mirna['alignment_list']:
                        nesteddict_gene_total[gene_id]['total_nr_interactions'] = nesteddict_gene_total[gene_id]['total_nr_interactions'] + 1
                        start = alignment['lnc_start']
                        end = alignment['lnc_end']
                        binding_site = start + "@@@" + end
                        interaction_bs = (gene_id, mirna_id, start, end)
                        if interaction_bs in found_bs:
                            #print("Existing interaction pair, existing binding site")
                            nesteddict[interaction_pair][binding_site] = nesteddict[interaction_pair][binding_site] + 1      
                        else:
                            #print("Existing interaction pair, new binding site")
                            found_bs.append(interaction_bs)
                            nesteddict[interaction_pair][binding_site] = 1
                            nesteddict_gene_total[gene_id]['total_unique_bs'] = nesteddict_gene_total[gene_id]['total_unique_bs'] + 1
                else:
                    #print("New interaction pair, new binding sites")
                    found_interaction_pair.append(interaction_pair)
                    for alignment in mirna['alignment_list']:
                        nesteddict_gene_total[gene_id]['total_nr_interactions'] = nesteddict_gene_total[gene_id]['total_nr_interactions'] + 1
                        start = alignment['lnc_start']
                        end = alignment['lnc_end']
                        binding_site = start + "@@@" + end
                        interaction_bs = (gene_id, mirna_id, start, end)
                        found_bs.append(interaction_bs)
                        nesteddict[interaction_pair][binding_site] = 1
                        nesteddict_gene_total[gene_id]['total_unique_bs'] = nesteddict_gene_total[gene_id]['total_unique_bs'] + 1

    for key, value in nesteddict.items():
        key_array = key.split("@@@")
        gene_id = key_array[0]
        mirna_id = key_array[1]
        unique_bs_per_interaction_pair = len(value.keys())
        total_interactions_per_pair = sum(value.values())

        unique_bs_per_gene_total = nesteddict_gene_total[gene_id]['total_unique_bs']
        total_interactions_per_gene = nesteddict_gene_total[gene_id]['total_nr_interactions']

        for bs, nr_interactions_current_bs in value.items():
            bs_array = bs.split("@@@")
            start = bs_array[0]
            stop = bs_array[1]
            output.write(gene_id + "\t" + mirna_id + "\t" + start + "\t" + stop + "\t" + str(nr_interactions_current_bs) + "\t"
                + str(unique_bs_per_interaction_pair) + "\t" + str(total_interactions_per_pair) + "\t" 
                + str(unique_bs_per_gene_total) + "\t" + str(total_interactions_per_gene) + "\n")


if __name__ == '__main__':
    cwd = sys.argv[1]
    fin = sys.argv[2]
    fout = sys.argv[3]

    f_out = os.path.join(cwd, fout)
    output = open(f_out, "w")

    nesteddict = defaultdict(dict)
    nesteddict_gene_total = defaultdict(dict)

    df = loadJson(fin)
    getInteractions(df)

# python3 getGeneMirnaInteractions2.py mm10/ mm10_interactions_allDBs_and_pc.json mm10_number_of_interactions.txt




