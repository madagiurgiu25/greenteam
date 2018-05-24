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

    fout=open(gtf + "out", "w")

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


if __name__ == "__main__":
    # asignUniqueKeys()

    fin = open("mapping_keys.txt","r").readlines()
    dict_mapping = {}
    for l in fin:
        arr = l.split("\t")
        if str(arr[0]).startswith("#") == False:
            dict_mapping[arr[1]] = arr[0]
    replaceKeyGTF("NONCODEv5_human_hg38_lncRNA.gtf",dict_mapping)
    replaceKeyGTF("NONCODEv5_mouse_mm10_lncRNA.gtf",dict_mapping)
    replaceKeyGTF("gencode.vM17.long_noncoding_RNAs.gtf",dict_mapping)
    replaceKeyGTF("gencode.v28.long_noncoding_RNAs.gtf",dict_mapping)
    # insertPKinJSON("hg38_long_noncoding_gencode.json",dict_mapping)
    # insertPKinJSON("hg38_long_noncoding_noncode.json",dict_mapping)
    # insertPKinJSON("mm10_long_noncoding_gencode.json",dict_mapping)
    # insertPKinJSON("mm10_long_noncoding_noncode.json",dict_mapping)

