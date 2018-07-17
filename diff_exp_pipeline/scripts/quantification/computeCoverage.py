__author__ = 'Mada'
# !/usr/bin/python3
# !C:\\Python34\\python.exe

########## find overlaps between genes
import json
import intervaltree
import re
import operator
from sympy import Interval, Union
import sys

# headers in GTF file
gtfdict = {'gene_id': "", 'transcript_id': "", 'transcript_type': "", 'transcript_name': "", 'gene_name': "",
           'gene_type': "", 'transcript_type': "", 'tag': "", 'exon_number': "", 'Parent': "", 'level': ""}

dict_words = {'geneID': "gene_id", 'transcriptID': "transcript_id", 'exonnumber': "exon_number"}

GENE_ID = 'gene_id'
TRANSCRIPT_ID = 'transcriptID'
EXON_NR = 'exonnumber'

EXON = 'exon'
GENE = 'gene'
TRANSCRIPT = 'transcript'
UTR = 'UTR'
CDS = 'CDS'

INTERVAL = 'interval'
COV = 'Coverage'
FPKM = 'FPKM'
TPM = 'TPM'
GENE_ID_STR = 'Gene ID'
GENE_NAME = 'Gene Name'
REF = 'Reference'
STRAND = 'Strand'
START = 'Start'
STOP = 'End'
REGION = 'Region'
dict_stringtie = {COV: "", FPKM: "", TPM: "", GENE_ID_STR: "", GENE_NAME: "", STRAND: "", START: "", STOP: ""}
RAW_COUNTS = 'Raw_counts'

# max number of gene overlapping partners (those parent genes will be ignored)
NMAX = 20


def union(data):
    """ Union of a list of intervals e.g. [(1,2),(3,4)] """
    intervals = [Interval(begin, end) for (begin, end) in data]
    u = Union(*intervals)
    # return u.measure
    return u
    # return [list(u.args[:2])] if isinstance(u, Interval) \
    #    else list(u.args)


def getGeneLength(gene):  # union exons

    list_exons = []

    for tr in gene["transcripts_list"]:
        for exon in tr["exon_list"]:
            list_exons.append((float(exon["start"]), float(exon["stop"])))

   
    union_list = union(list_exons)
    print(gene["gene_id"])
    
    count_len = union_list.measure + len(list(union_list.args))
    
    return int(count_len)


def computePerSample(jsonGTF, geneabundanceFile):
    count = 0
    arr_header = []

    dict_abundance = {}

    # load existing anno file
    try:
        with open(jsonGTF, 'r') as f:
            jsonGTF_genes = json.load(f)
    except FileNotFoundError:
        print("your file does not exist")

    with open(geneabundanceFile, 'r', encoding='utf-8') as infile:
        for line in infile:
            if count == 0:
                # read header
                arr_header = line.strip().replace("\n", "").split("\t")
                count = 1
            else:
                arr = line.strip().replace("\n", "").split("\t")
                dict_aux = {}
                for i, val in enumerate(arr):
                    dict_aux[arr_header[i]] = val

                if arr[0] not in dict_abundance:
                    dict_abundance[arr[0]] = {}
                    dict_abundance[arr[0]][INTERVAL] = []
                    dict_abundance[arr[0]][FPKM] = 0
                    dict_abundance[arr[0]][TPM] = 0
                    dict_abundance[arr[0]][COV] = 0
                dict_abundance[arr[0]][FPKM] += float(dict_aux[FPKM])
                dict_abundance[arr[0]][TPM] += float(dict_aux[TPM])
                dict_abundance[arr[0]][COV] += float(dict_aux[COV])
                dict_abundance[arr[0]][STRAND] = dict_aux[STRAND]
                dict_abundance[arr[0]][GENE_NAME] = dict_aux[GENE_NAME]
                dict_abundance[arr[0]][INTERVAL].append(float(arr[4]))
                dict_abundance[arr[0]][INTERVAL].append(float(arr[5]))
                dict_abundance[arr[0]][REF] = dict_aux[REF]
                dict_abundance[arr[0]][START] = dict_aux[START]
                dict_abundance[arr[0]][STOP] = dict_aux[STOP]

        for gene in jsonGTF_genes:
            if gene[GENE_ID] in dict_abundance:
                dict_abundance[gene[GENE_ID]][RAW_COUNTS] = gene["length"] * dict_abundance[gene[GENE_ID]][COV]
                dict_abundance[gene[GENE_ID]][START] = min(dict_abundance[gene[GENE_ID]][INTERVAL])
                dict_abundance[gene[GENE_ID]][STOP] = max(dict_abundance[gene[GENE_ID]][INTERVAL])

        with open(geneabundanceFile + "_rawcounts", 'w') as outfile:
            outfile.write("\t".join(str(item) for item in arr_header) + "\t" + RAW_COUNTS + "\n")
            for gene in dict_abundance:
                if RAW_COUNTS not in dict_abundance[gene]:
                    outfile.write(gene + "\t" +
                                  dict_abundance[gene][GENE_NAME] + "\t" +
                                  dict_abundance[gene][REF] + "\t" +
                                  dict_abundance[gene][STRAND] + "\t" +
                                  str(dict_abundance[gene][START]) + "\t" +
                                  str(dict_abundance[gene][STOP]) + "\t" +
                                  str(dict_abundance[gene][COV]) + "\t" +
                                  str(dict_abundance[gene][FPKM]) + "\t" +
                                  str(dict_abundance[gene][TPM]) + "\t" +
                                  str(0) + "\n")
                else:
                    outfile.write(gene + "\t" +
                                  dict_abundance[gene][GENE_NAME] + "\t" +
                                  dict_abundance[gene][REF] + "\t" +
                                  dict_abundance[gene][STRAND] + "\t" +
                                  str(dict_abundance[gene][START]) + "\t" +
                                  str(dict_abundance[gene][STOP]) + "\t" +
                                  str(dict_abundance[gene][COV]) + "\t" +
                                  str(dict_abundance[gene][FPKM]) + "\t" +
                                  str(dict_abundance[gene][TPM]) + "\t" +
                                  str(dict_abundance[gene][RAW_COUNTS]) + "\n")

def compute(jsonGTF, geneabundancelist):  # json gtf in house + geneabundance file from stringtie

    count = 0
    arr_header = []

    dict_abundance = {}

    # load existing anno file
    try:
        with open(jsonGTF, 'r') as f:
            jsonGTF_genes = json.load(f)
    except FileNotFoundError:
        print("your file does not exist")

    for file in geneabundancelist:
        print("process file: " + file)
        # parse gene_abundance entry to dict
        count = 0
        with open(file, 'r', encoding='utf-8') as infile:
            for line in infile:
                if count == 0:
                    # read header
                    arr_header = line.strip().replace("\n", "").split("\t")
                    count = 1
                else:
                    arr = line.strip().replace("\n", "").split("\t")
                    dict_aux = {}
                    for i, val in enumerate(arr):
                        dict_aux[arr_header[i]] = val

                    if arr[0] not in dict_abundance:
                        dict_abundance[arr[0]] = {}
                        dict_abundance[arr[0]][INTERVAL] = []
                        dict_abundance[arr[0]][FPKM] = 0
                        dict_abundance[arr[0]][TPM] = 0
                        dict_abundance[arr[0]][COV] = 0
                    dict_abundance[arr[0]][FPKM] += float(dict_aux[FPKM])
                    dict_abundance[arr[0]][TPM] += float(dict_aux[TPM])
                    dict_abundance[arr[0]][COV] += float(dict_aux[COV])
                    dict_abundance[arr[0]][STRAND] = dict_aux[STRAND]
                    dict_abundance[arr[0]][GENE_NAME] = dict_aux[GENE_NAME]
                    dict_abundance[arr[0]][INTERVAL].append(float(arr[4]))
                    dict_abundance[arr[0]][INTERVAL].append(float(arr[5]))
                    dict_abundance[arr[0]][REF] = dict_aux[REF]
                    dict_abundance[arr[0]][START] = dict_aux[START]
                    dict_abundance[arr[0]][STOP] = dict_aux[STOP]

        for gene in jsonGTF_genes:
            if gene[GENE_ID] in dict_abundance:
                dict_abundance[gene[GENE_ID]][RAW_COUNTS] = gene["length"] * dict_abundance[gene[GENE_ID]][COV]
                dict_abundance[gene[GENE_ID]][START] = min(dict_abundance[gene[GENE_ID]][INTERVAL])
                dict_abundance[gene[GENE_ID]][STOP] = max(dict_abundance[gene[GENE_ID]][INTERVAL])

        with open(file + "_rawcounts", 'w') as outfile:
            outfile.write("\t".join(str(item) for item in arr_header) + "\t" + RAW_COUNTS + "\n")
            for gene in dict_abundance:
                if RAW_COUNTS not in dict_abundance[gene]:
                    outfile.write(gene + "\t" +
                                  dict_abundance[gene][GENE_NAME] + "\t" +
                                  dict_abundance[gene][REF] + "\t" +
                                  dict_abundance[gene][STRAND] + "\t" +
                                  str(dict_abundance[gene][START]) + "\t" +
                                  str(dict_abundance[gene][STOP]) + "\t" +
                                  str(dict_abundance[gene][COV]) + "\t" +
                                  str(dict_abundance[gene][FPKM]) + "\t" +
                                  str(dict_abundance[gene][TPM]) + "\t" +
                                  str(0) + "\n")
                else:
                    outfile.write(gene + "\t" +
                                  dict_abundance[gene][GENE_NAME] + "\t" +
                                  dict_abundance[gene][REF] + "\t" +
                                  dict_abundance[gene][STRAND] + "\t" +
                                  str(dict_abundance[gene][START]) + "\t" +
                                      str(dict_abundance[gene][STOP]) + "\t" +
                                  str(dict_abundance[gene][COV]) + "\t" +
                                  str(dict_abundance[gene][FPKM]) + "\t" +
                                  str(dict_abundance[gene][TPM]) + "\t" +
                                  str(dict_abundance[gene][RAW_COUNTS]) + "\n")


def computeGeneLength(jsonGTF):
    # load existing anno file
    try:
        with open(jsonGTF, 'r') as f:
            jsonGTF_genes = json.load(f)
    except FileNotFoundError:
        print("your file does not exist")

    for gene in jsonGTF_genes:
        gene["length"] = getGeneLength(gene)

    with open("hg38_primary_assembly_and_lncRNA_withgenelength.json", 'w') as outfile:
        json.dump(jsonGTF_genes, outfile, indent=4)


if __name__ == "__main__":

    #geneabundancefile = ["M103A1/M103A1_geneabundance.tab", "M103A2/M103A2_geneabundance.tab",
    #                     "M103A3/M103A3_geneabundance.tab", "M103A4/M103A4_geneabundance.tab",
    #                     "M103A5/M103A5_geneabundance.tab", "M103B1/M103B1_geneabundance.tab",
    #                     "M103B2/M103B2_geneabundance.tab", "M103B3/M103B3_geneabundance.tab",
    #                     "M103B4/M103B4_geneabundance.tab", "M103B5/M103B5_geneabundance.tab",
    #                     "Mlet7A1/Mlet7A1_geneabundance.tab", "Mlet7A2/Mlet7A2_geneabundance.tab",
    #                     "Mlet7A3/Mlet7A3_geneabundance.tab", "Mlet7A4/Mlet7A4_geneabundance.tab",
    #                     "Mlet7A7/Mlet7A7_geneabundance.tab", "Mlet7A8/Mlet7A8_geneabundance.tab",
    #                     "Mlet7A9/Mlet7A9_geneabundance.tab", "Mlet7B1/Mlet7B1_geneabundance.tab",
    #                     "Mlet7B2/Mlet7B2_geneabundance.tab", "Mlet7B3/Mlet7B3_geneabundance.tab",
    #                     "Mlet7B4/Mlet7B4_geneabundance.tab", "Mlet7B7/Mlet7B7_geneabundance.tab",
    #                     "Mlet7B8/Mlet7B8_geneabundance.tab", "Mlet7B9/Mlet7B9_geneabundance.tab",
    #                     ]

    #compute("mm10_primary_assembly_and_lncRNA_withgenelength.json", geneabundancefile)

    if len(sys.argv) == 4:
        jsonGTF = sys.argv[1]
        geneabundancefile = sys.argv[2]
        outfile = sys.argv[3]
        compute(jsonGTF, geneabundancefile, outfile)
   
    # # add gene length to dictionary
    elif len(sys.argv) == 2:
        computeGeneLength(sys.argv[1])
