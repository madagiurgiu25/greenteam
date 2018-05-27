#!/usr/bin/python3
#!C:\\Python34\\python.exe

import json
count = 0
ALIAS = 'alias'
SPECIES = 'species'
TYPE = 'type'

def addKeys(file, type, type_long, species, dict ):
    global count
    f1 = open(file, 'r').readlines()
    for l in f1:
        val = l.strip()
        if len(val) > 0 and val.startswith("#") == False: # jump over empty entries'
            key = "LNC_" + type + "_" + species + "_" + str(count).zfill(8)
            dict[key] = {}
            dict[key][ALIAS] = val
            dict[key][SPECIES] = species
            dict[key][TYPE] = type_long
            count = count + 1

    return dict

def asignUniqueKeys():
    dict_keys={}

    dict_keys = addKeys("genes_hg38.txt", "GE", "GENE", "hg38", dict_keys)
    dict_keys = addKeys("transcripts_hg38.txt", "TR", "TRANSCRIPT", "hg38", dict_keys)
    dict_keys = addKeys("genes_mm10.txt", "GE","GENE", "mm10", dict_keys)
    dict_keys = addKeys("transcripts_mm10.txt", "TR", "TRANSCRIPT", "mm10", dict_keys)

    with open("mapping_keys.json", 'w') as outfile:
        json.dump(dict_keys, outfile, indent=4)
    outfile.close()

    fout = open("mapping_keys.txt","w")
    fout.write("#mingleID\talias\tspecies\ttype\n")
    for key in dict_keys:
        fout.write(key + "\t" + dict_keys[key][ALIAS] + "\t" + dict_keys[key][SPECIES] + "\t" + dict_keys[key][TYPE] + "\n")

    fout.close()

def insertPKinJSON(file):

    pass

connector = " "
def replaceKeyGTF(gtf,mapObj):

    fout=open(gtf + "_newkeys", "w")

    try:
        with open(gtf, "r") as f:
            for line in f:
                if str(line).startswith("#"):
                    fout.write(line)
                else:
                    info_head = line.strip().replace("\n","").split("\t")[0:8]
                    print(info_head)
                    info_row = line.strip().replace("\n","").split("\t")[8]
                    info_arr = info_row.replace("\"","").split(";")
                    info_tail = []
                    for item in info_arr:
                        if len(item.strip())>0:
                            (key,val) = item.strip().split(" ")
                            if key in ["gene_id","transcript_id"]:
                                if val in mapObj:
                                    info_tail.append(key + connector + "\"" + mapObj[val] + "\"" )
                                    info_tail.append(key + "_alias" + connector + "\"" + val + "\"")
                                else:
                                    info_tail.append(key + connector + "\"" + val + "\"")
                            else:
                                info_tail.append(key + connector + "\"" + val + "\"")
                    string_tail = "; ".join(info_tail)
                    info_head.append(string_tail)
                    fout.write("\t".join(info_head) + "\n")
        fout.close()

    except FileNotFoundError:
        print("your file does not exist")

def extractFPKM_TPM():
    pass

def addGeneEntrie(gtf):

    list_trascript_start_stop = []
    current_gene_id = ""
    current_strand = ""
    current_chr = ""
    try:
        fout=open(gtf + "_withgenes", "w")
        with open(gtf, "r") as f:
            for line in f:

                if line.startswith("#"):
                    fout.write(line)
                    continue
                else:
                    arr = line.strip().replace("\n","").split("\t")
                    if arr[2] != 'transcript':
                        fout.write(line)
                        continue
                    # check if transcript on same gene and save start stop
                    info_row = arr[8]
                    start = int(arr[3])
                    stop = int(arr[4])
                    info_arr = info_row.replace("\"","").split(";")
                    for item in info_arr:
                        if len(item.strip())>0:
                            (key,val) = item.strip().split(" ")
                            if key  == "gene_id":
                                if val == current_gene_id: # means that transcript its on the same gene
                                    list_trascript_start_stop.append(start)
                                    list_trascript_start_stop.append(stop)
                                    current_strand = arr[6]
                                    current_chr = arr[0]
                                else:
                                    if len(current_gene_id) > 0: # do i have already a gene in memory
                                        max_stop = max(list_trascript_start_stop)
                                        min_start = min(list_trascript_start_stop)

                                        fout.write(("{0}\tCufflinks\tgene\t{1}\t{2}\t0\t{3}\t.\tgene_id \"{4}\";\n").format(current_chr,min_start,max_stop,current_strand,current_gene_id))

                                    list_trascript_start_stop = []
                                    list_trascript_start_stop.append(start)
                                    list_trascript_start_stop.append(stop)
                                    current_gene_id = val
                                    current_strand = arr[6]
                                    current_chr = arr[0]
                                break
                    fout.write(line)

            if len(current_gene_id) > 0:
                max_stop = max(list_trascript_start_stop)
                min_start = min(list_trascript_start_stop)

                fout.write(("{0}\tCufflinks\tgene\t{1}\t{2}\t0\t{3}\t.\tgene_id \"{4}\";\n").format(current_chr,min_start,max_stop,current_strand,current_gene_id))
            fout.close()

    except FileNotFoundError:
        print("your file does not exist")

if __name__ == "__main__":

    ################ Assgin unique keys #####################
    #asignUniqueKeys()

    ################ Add gene entries #####################
    # addGeneEntrie("NONCODEv5_mouse_mm10_lncRNA.gtf")
    # addGeneEntrie("NONCODEv5_human_hg38_lncRNA.gtf")
    addGeneEntrie("lncrnadb_mm10.gtf")


    ################ Mapping keyes #####################
    fin = open("mapping_keys.txt","r").readlines()
    dict_mapping = {}
    for l in fin:
        arr = l.split("\t")
        if str(arr[0]).startswith("#") == False:
            dict_mapping[arr[1]] = arr[0]
    replaceKeyGTF("lncrnadb_mm10.gtf_withgenes",dict_mapping)
    #replaceKeyGTF("NONCODEv5_mouse_mm10_lncRNA.gtf_withgenes",dict_mapping)

    #replaceKeyGTF("gencode.vM17.long_noncoding_RNAs.gtf",dict_mapping)
    #replaceKeyGTF("gencode.v28.long_noncoding_RNAs.gtf",dict_mapping)

