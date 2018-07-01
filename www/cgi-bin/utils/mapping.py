#!/usr/bin/python3
#!C:\\Python34\\python.exe

import json
import re
count = 0
ALIAS = 'alias'
SPECIES = 'species'
TYPE = 'type'
CHR = 'chr'
SOURCE = 'source'
ID = 'id'
PARENT_GENE = "parent_gene_id"
PARENT_GENE_ALIAS = "parent_gene_alias"
GENE = "GENE"
TRANSCRIPT = 'TRANSCRIPT'
COUNT = "COUNT"

GENCODE='gencode'
LNCIPEDIA='lncipedia'
NONCODE='noncode'
LNCRNADB='lncrnadb'

dict_positions={GENCODE:0,LNCIPEDIA:1,NONCODE:2,LNCRNADB:3}

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
def replaceKeyGTF(gtf,type):

    try:
        fout=open(gtf + "_newkeys", "w")

        fin = open("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/mapping/lncRNA/mapping_keys.txt","r")
        mapObj = {}
        mapObj[TRANSCRIPT] = {}
        mapObj[GENE] = {}
        for l in fin:
            if l.startswith("#") == False:
                arr = l.strip().replace("\n","").split("\t")
                if arr[2] == type:
                    mapObj[arr[3]][arr[1]] = arr[0]
        fin.close()
        # fmap = open("../mapping_keys.txt", "a")
        # count_tr = 733858
        #
        # with open('../mapping_keys.json') as f:
        #     mapObjNew = json.load(f)

        with open(gtf, "r") as f:
            for line in f:
                if str(line).startswith("#"):
                    fout.write(line)
                else:
                    info_head = line.strip().replace("\n","").split("\t")[0:8]
                    # print(info_head)
                    info_row = line.strip().replace("\n","").split("\t")[8]
                    info_arr = info_row.replace("\"","").split(";")
                    info_tail = []
                    for item in info_arr:
                        if len(item.strip())>0:
                            splitvals = item.strip().split(" ")
                            if len(splitvals) == 2:
                                (key,val) = splitvals
                                if key in ["gene_id","transcript_id"]:
                                    # if key == "gene_id":
                                    #     type = "GE"
                                    #     type_long = "GENE"
                                    # else:
                                    #     type = "TR"
                                    #     type_long = "TRANSCRIPT"
                                    # if val not in mapObj: # assign extra keys
                                    #     key_map = "LNC_" + type + "_" + "hg38" + "_" + str(count_tr).zfill(8)
                                    #     mapObjNew[key_map] = {}
                                    #     mapObjNew[key_map][ALIAS] = val
                                    #     mapObjNew[key_map][SPECIES] = "hg38"
                                    #     mapObjNew[key_map][TYPE] = type_long
                                    #     count_tr = count_tr + 1
                                    #     print(key_map + " " + val + " " + type + " " + type_long)
                                    #     fmap.write(key_map   + "\t" + mapObjNew[key_map][ALIAS] + "\t" + mapObjNew[key_map][SPECIES] + "\t" + mapObjNew[key_map][TYPE] + "\n")
                                    #     mapObj[val] = key_map

                                    # find LNC key
                                    if key == 'transcript_id' and val in mapObj[TRANSCRIPT]:
                                        info_tail.append(key + connector + "\"" + mapObj[TRANSCRIPT][val] + "\"" )
                                        info_tail.append(key + "_alias" + connector + "\"" + val + "\"")
                                    elif key == 'gene_id' and val in mapObj[GENE]:
                                        info_tail.append(key + connector + "\"" + mapObj[GENE][val] + "\"" )
                                        info_tail.append(key + "_alias" + connector + "\"" + val + "\"")
                                    else:
                                        print(key + " " + val + " not found!!")
                                        info_tail.append(key + connector + "\"" + val + "\"")
                                else:
                                    info_tail.append(key + connector + "\"" + val + "\"")
                            else:
                                print("error with space " + item)
                    string_tail = "; ".join(info_tail)
                    info_head.append(string_tail)
                    fout.write("\t".join(info_head) + "\n")
        fout.close()
        # fmap.close()
        # with open("../mapping_keys.json", 'w') as outfile:
        #     json.dump(mapObjNew, outfile, indent=4)
        # outfile.close()
    except FileNotFoundError:
        print("your file does not exist")

def extractFPKM_TPM():
    pass

def addGeneEntrie(gtf):

    dict_genes = {}
    list_trascript_start_stop = "list_intervals"
    current_gene_id = ""
    current_strand = ""
    current_chr = ""
    chr="chr"
    strand="strand"

    try:
        fout=open(gtf + "_withgenes", "w")
        with open(gtf, "r") as f:
            for line in f:

                if line.startswith("#"):
                    fout.write(line)
                    continue
                else:
                    arr = line.strip().replace("\n","").split("\t")
                    if arr[2] != TRANSCRIPT:
                        fout.write(line)
                        continue

                    fout.write(line)
                    # check if transcript on same gene and save start stop
                    info_row = arr[8]
                    start = int(arr[3])
                    stop = int(arr[4])
                    current_strand = arr[6]
                    current_chr = arr[0]
                    info_arr = info_row.replace("\"","").split(";")
                    for item in info_arr:
                        if len(item.strip())>0:
                            (key,val) = item.strip().split(" ")
                            if key  == "gene_id":
                                if val in dict_genes:
                                    dict_genes[val][list_trascript_start_stop].append(start)
                                    dict_genes[val][list_trascript_start_stop].append(stop)
                                else:
                                    dict_genes[val] = {}
                                    dict_genes[val][chr]=current_chr
                                    dict_genes[val][strand]=current_strand
                                    dict_genes[val][list_trascript_start_stop] = []
                                    dict_genes[val][list_trascript_start_stop].append(start)
                                    dict_genes[val][list_trascript_start_stop].append(stop)

                                break

            for k in dict_genes:
                current_chr=dict_genes[k][chr]
                current_strand=dict_genes[k][strand]
                current_gene_id=k
                max_stop = max(dict_genes[k][list_trascript_start_stop])
                min_start = min(dict_genes[k][list_trascript_start_stop])

                fout.write(("{0}\tCufflinks\tgene\t{1}\t{2}\t0\t{3}\t.\tgene_id \"{4}\";\n").format(current_chr,min_start,max_stop,current_strand,current_gene_id))
            fout.close()

    except FileNotFoundError:
        print("your file does not exist")


def addGeneTranscriptEntry(gtf):

    dict_genes = {}
    list_trascript_start_stop = "list_intervals"
    list_exons_start_stop = "list_exons"
    transcript_list = "transcript_list"
    current_gene_id = ""
    current_strand = ""
    current_chr = ""
    chr="chr"
    strand="strand"

    try:
        fout=open(gtf + "_withgenes", "w")
        with open(gtf, "r") as f:
            for line in f:

                if line.startswith("#"):
                    fout.write(line)
                    continue
                else:
                    arr = line.strip().replace("\n","").split("\t")

                    fout.write(line)
                    # check if transcript on same gene and save start stop
                    info_row = arr[8]
                    start = int(arr[3])
                    stop = int(arr[4])
                    current_strand = arr[6]
                    current_chr = arr[0]
                    info_arr = info_row.replace("\"","").split(";")

                    geneid = re.search("gene_id \"(.*?)\" ;", arr[8], flags=0)
                    if geneid:
                        geneid = geneid.group(1)

                    transcriptid = re.search("transcript_id \"(.*?)\" ;", arr[8], flags=0)
                    if transcriptid:
                        transcriptid = transcriptid.group(1)

                    print(str(geneid) + "\t" + str(transcriptid))

                    if geneid in dict_genes:
                        dict_genes[geneid][list_trascript_start_stop].append(start)
                        dict_genes[geneid][list_trascript_start_stop].append(stop)
                    else:
                        dict_genes[geneid] = {}
                        dict_genes[geneid][chr]=current_chr
                        dict_genes[geneid][strand]=current_strand
                        dict_genes[geneid][list_trascript_start_stop] = []
                        dict_genes[geneid][list_trascript_start_stop].append(start)
                        dict_genes[geneid][list_trascript_start_stop].append(stop)
                        dict_genes[geneid][transcript_list] = {}

                    if transcriptid in dict_genes[geneid][transcript_list]:
                        dict_genes[geneid][transcript_list][transcriptid][list_exons_start_stop].append(start)
                        dict_genes[geneid][transcript_list][transcriptid][list_exons_start_stop].append(stop)
                    else:
                        dict_genes[geneid][transcript_list][transcriptid] = {}
                        dict_genes[geneid][transcript_list][transcriptid][chr] = current_chr
                        dict_genes[geneid][transcript_list][transcriptid][strand] = current_strand
                        dict_genes[geneid][transcript_list][transcriptid][list_exons_start_stop] = []
                        dict_genes[geneid][transcript_list][transcriptid][list_exons_start_stop].append(start)
                        dict_genes[geneid][transcript_list][transcriptid][list_exons_start_stop].append(stop)

            for k in dict_genes:
                current_chr=dict_genes[k][chr]
                current_strand=dict_genes[k][strand]
                current_gene_id=k
                max_stop = max(dict_genes[k][list_trascript_start_stop])
                min_start = min(dict_genes[k][list_trascript_start_stop])

                fout.write(("{0}\tCufflinks\tgene\t{1}\t{2}\t0\t{3}\t.\tgene_id \"{4}\";\n").format(current_chr,min_start,max_stop,current_strand,current_gene_id))
                for t in dict_genes[k][transcript_list]:
                    current_chr=dict_genes[k][transcript_list][t][chr]
                    current_strand=dict_genes[k][transcript_list][t][strand]
                    current_transcript_id=k
                    max_stop = max(dict_genes[k][transcript_list][t][list_exons_start_stop])
                    min_start = min(dict_genes[k][transcript_list][t][list_exons_start_stop])
                    fout.write(("{0}\tCufflinks\ttranscript\t{1}\t{2}\t0\t{3}\t.\tgene_id \"{4}\"; transcript_id \"{5}\";\n").format(current_chr,min_start,max_stop,current_strand,current_gene_id,current_transcript_id))
            fout.close()

    except FileNotFoundError:
        print("your file does not exist")

# transform a json structure in a dict like "transcript_id":"gene_id"
# json and txt
def transcriptAndGeneDict(jsonFile,source,outfile):
    dict = {}
    try:
        fout=open(outfile, "a")
        with open(jsonFile) as f:
            mapObjNew = json.load(f)
            for ent in mapObjNew:
                for t in ent["transcripts_list"]:
                    dict[t["transcript_id"]] = ent["gene_id"]
                    fout.write(ent["gene_id"] + "\t" + t["transcript_id"] + "\t" + source + "\n")
        fout.close()
        return dict
    except FileNotFoundError:
        print("your file does not exist")

def mergeNewKeysInOverlap(dict_mapping_NAME2ID,dict_TR_GE,dict_overlap,species,dict_mapping_ID2NAME):

    f = open("overlap_" + species + "_newkeys_complete.txt","w")
    f.write("#lnc_overlap_key\tgencode_tr\tlncipedia_tr\tnoncode_tr\tlncrnadb_tr\tgencode_ge\tlncipedia_ge\tnoncode_ge\tlncrnadb_ge\n")
    for ent in dict_overlap:
        for key in ent: # overlap information
            arr = ["","","",""]
            arr_keys = ["","","",""]
            for elem in ent[key]:
                elem[ID] = dict_mapping_NAME2ID[TRANSCRIPT][elem[ALIAS]]
                elem[PARENT_GENE] = dict_TR_GE[elem[SOURCE]][elem[ID]]
                elem[PARENT_GENE_ALIAS] = dict_mapping_ID2NAME[GENE][elem[PARENT_GENE]]

                arr[dict_positions[elem[SOURCE]]] = elem[ID]
                arr_keys[dict_positions[elem[SOURCE]]] = elem[PARENT_GENE]
            f.write(key + "\t" + arr[0] + "\t" + arr[1] + "\t" + arr[2] + "\t" + arr[3] + "\t" + arr_keys[0] + "\t" + arr_keys[1] + "\t" + arr_keys[2] + "\t" + arr_keys[3] + "\n")
    f.close()
    return dict_overlap

def build_nonredundantset(species):

    bigjsonGeneTranscript = {}
    bigjsonGeneTranscript[LNCIPEDIA] = transcriptAndGeneDict(species + "_lncipedia_newkeys.json","lncipedia",species + "_mapGeneTranscripts.txt")
    bigjsonGeneTranscript[NONCODE] = transcriptAndGeneDict(species + "_noncode_newkeys.json","noncode",species + "_mapGeneTranscripts.txt")
    bigjsonGeneTranscript[LNCRNADB] =transcriptAndGeneDict(species + "_lncrnadb_newkeys.json","lncrnadb",species + "_mapGeneTranscripts.txt")
    bigjsonGeneTranscript[GENCODE] = transcriptAndGeneDict(species + "_gencode_long_noncoding_newkeys.json","gencode",species + "_mapGeneTranscripts.txt")

    with open(species + "_mapGeneTranscripts.json", 'w') as outfile:
        json.dump(bigjsonGeneTranscript, outfile, indent=4)

    ################## complete the overlapping dict
    fin = open("../mapping_keys.txt","r").readlines()
    dict_mapping_NAME2ID = {} # gene name to LNC id
    dict_mapping_NAME2ID[TRANSCRIPT] = {}
    dict_mapping_NAME2ID[GENE] = {}

    dict_mapping_ID2NAME = {} # gene name to LNC id
    dict_mapping_ID2NAME[TRANSCRIPT] = {}
    dict_mapping_ID2NAME[GENE] = {}

    for l in fin:
        arr = l.replace("\n","").split("\t")
        if str(arr[0]).startswith("#") == False and arr[2] == species:
            dict_mapping_NAME2ID[arr[3]][arr[1]] = arr[0]
            dict_mapping_ID2NAME[arr[3]][arr[0]] = arr[1]

    with open(species + "_mapGeneTranscripts.json", 'r') as f:
        dict_TR_GE = json.load(f)

    with open("overlap_" + species + ".json", 'r') as f:
        dict_overlap = json.load(f)

    dict_overlap_complete = mergeNewKeysInOverlap(dict_mapping_NAME2ID,dict_TR_GE,dict_overlap,species,dict_mapping_ID2NAME)
    with open("overlap_" + species + "_newkeys_complete.json", 'w') as outfile:
        json.dump(dict_overlap_complete, outfile, indent=4)


def insertKeys(mapping,jsonmapping,newkeylist,type,species,typeshort):

    dict_new_keys = {}
    fopen=open(newkeylist,"r").readlines()
    for f in fopen:
        k = f.strip().replace("\n","")
        dict_new_keys[k] = {}
        dict_new_keys[k][SPECIES] = species
        dict_new_keys[k][TYPE] = type

    last = None
    fopen=open(mapping,"r").readlines()
    for f in fopen:
        line = f.strip().replace("\n","").split("\t")
        last = line[0] # letzte elem

        # remove all entries from dict which are already mapped
        if line[1] in dict_new_keys and line[2] == species and line[3] == type:
            del dict_new_keys[line[1]]

    count = int(last.split("_")[3]) + 1
    print("Old count:" + str(count))
    if len(dict_new_keys) > 0:
        fout=open(mapping,"a")

        with open(jsonmapping, 'r') as infile:
            dict_mapped_keys = json.load(infile)

        for key in dict_new_keys:
            key_lnc = "LNC_" + typeshort + "_" + species + "_" + str(count).zfill(8)
            count = count + 1
            fout.write(key_lnc + "\t" + key + "\t" + dict_new_keys[key][SPECIES] + "\t" + dict_new_keys[key][TYPE] + "\n")
            dict_mapped_keys[key_lnc] = {}
            dict_mapped_keys[key_lnc][ALIAS] = key
            dict_mapped_keys[key_lnc][SPECIES] = dict_new_keys[key][SPECIES]
            dict_mapped_keys[key_lnc][TYPE] = dict_new_keys[key][TYPE]
            print(key_lnc + "\t" + key + "\t" + dict_new_keys[key][SPECIES] + "\t" + dict_new_keys[key][TYPE] + "\n")
        fout.close()

        with open(jsonmapping, 'w') as outfile:
            json.dump(dict_mapped_keys, outfile, indent=4)

def removeRedundancy(species):

    tr_included = {}
    tr_excluded = {}
    gene_included = {}
    gene_excluded = {}
    gene_unclear = {}
    COUNT_EX = 'count_ex'
    COUNT_IN = 'count_in'
    count_conflict = 0


    with open("overlap_" + species + "_newkeys_complete.json", 'r') as f:
        dict_overlap_complete = json.load(f)

    with open(species + "_mapGeneTranscripts.json", 'r') as f:
        dict_TR_GE = json.load(f)

    dict_GE_TR = {}
    for source in dict_TR_GE:
        if source not in dict_GE_TR:
            dict_GE_TR[source] = {}
        for m in dict_TR_GE[source]:
            key_gene = dict_TR_GE[source][m]
            if key_gene  not in dict_GE_TR[source]: # gene not in dictionary
                dict_GE_TR[source][key_gene] = {}
                dict_GE_TR[source][key_gene][TRANSCRIPT] = []
                dict_GE_TR[source][key_gene][COUNT] = 0
            dict_GE_TR[source][key_gene][TRANSCRIPT].append(m)
            dict_GE_TR[source][key_gene][COUNT] += 1

    for item in dict_overlap_complete:
        for lnc_consensus in item:
            # if len(item[lnc_consensus]) == 1:
            #     for overlap_item in item[lnc_consensus]:
            #         tr_included[overlap_item[ID]] = dict_TR_GE[overlap_item[SOURCE]][overlap_item[ID]]
            #         gene_included[overlap_item[PARENT_GENE]] = dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]]
            #         break
            #     continue

            if len(item[lnc_consensus]) > 1: # least one overlap
                min_transcripts = 100000
                max_transcripts = 0
                include_tr = ""

                # find the gene in the overlaps which has the min number of transcripts
                for overlap_item in item[lnc_consensus]:
                    if overlap_item[PARENT_GENE] in gene_included:
                        include_tr = overlap_item[ID]
                        break
                    if overlap_item[PARENT_GENE] not in gene_excluded and dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][COUNT] < min_transcripts: # found a gene already included
                        include_tr = overlap_item[ID]
                        min_transcripts = dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][COUNT]

                    # if overlap_item[PARENT_GENE] not in gene_excluded and dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][COUNT] > max_transcripts: # found a gene already included
                    #     include_tr = overlap_item[ID]
                    #     max_transcripts = dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][COUNT]

                if include_tr != "":
                    # sort the transcripts and their genes in the 2 lists (included/excluded)
                    for overlap_item in item[lnc_consensus]:
                        if include_tr != overlap_item[ID]:
                            if overlap_item[PARENT_GENE] not in gene_included: # gene is not on included list
                                if overlap_item[PARENT_GENE] not in gene_excluded: # gene is excluded
                                    # print("exclude gene: " + overlap_item[PARENT_GENE])
                                    gene_excluded[overlap_item[PARENT_GENE]] = ''
                                for tr in dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][TRANSCRIPT]: # remove all transcripts
                                    if tr not in tr_excluded:
                                        tr_excluded[tr] = ""
                                    if tr == overlap_item[ID]:
                                        tr_excluded[tr] = include_tr
                                        # print("and its transcripts :" + tr)
                                # else:
                                    # print("already excluded gene: " + overlap_item[PARENT_GENE])
                            else:
                                # print("gene already included")
                                # gene is on the included list
                                tr_included[overlap_item[ID]] = ''
                        else:
                            if overlap_item[PARENT_GENE] not in gene_included: # gene is not on unclear state
                                # print("include gene: " + overlap_item[PARENT_GENE])
                                gene_included[overlap_item[PARENT_GENE]] = ''
                            for tr in dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][TRANSCRIPT]: # remove all transcripts
                                if tr not in tr_included:
                                    tr_included[tr] = ''
                                # print("and its transcripts :" + tr)
                            # else:
                                # print("already included gene: " + overlap_item[PARENT_GENE])
                else: # all genes are on excluded
                    # solve conflict
                    max_empty_transcripts = 0
                    include_tr = ""
                    count_conflict +=1
                    # print(item[lnc_consensus])
                    for overlap_item in item[lnc_consensus]: # find the gene with the most not covered transcripts
                        max_local = 0
                        if tr_excluded[overlap_item[ID]] == "":
                            for tr in dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][TRANSCRIPT]:
                                if tr_excluded[tr] == "":
                                    max_local += 1
                            if max_local > max_empty_transcripts:
                                include_tr = overlap_item[ID]
                                max_empty_transcripts = max_local
                    if include_tr != "": # remove from excluded add to included the gene with the highest nr of transcripts uncovered
                        for overlap_item in item[lnc_consensus]:
                            if include_tr == overlap_item[ID]:
                                gene_included[overlap_item[PARENT_GENE]] = ''
                                del gene_excluded[overlap_item[PARENT_GENE]]
                                for tr in dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][TRANSCRIPT]: # remove all transcripts
                                    if tr not in tr_included:
                                        tr_included[tr] = ''
                                    del tr_excluded[tr]

    # solve single transcript...not overlapping
    # if transcript sitting on an included gene everything fine
    # if parent gene nor in excluded or included everything fine
    # if parent gene on excluded check the transcripts...if there are transcript not excluded by another transcript put the gene on unclear
    for item in dict_overlap_complete:
        for lnc_consensus in item:
            if len(item[lnc_consensus]) == 1: # transcript does not overlap anything
                for overlap_item in item[lnc_consensus]:
                    if overlap_item[PARENT_GENE] not in gene_included and overlap_item[PARENT_GENE] not in gene_excluded:
                        # print("found single transcript: include " + overlap_item[PARENT_GENE] + " " + overlap_item[ID])
                        gene_included[overlap_item[PARENT_GENE]] = ''
                        for tr in dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][TRANSCRIPT]: # remove all transcripts
                            if tr not in tr_included:
                                tr_included[tr] = ''
                    elif overlap_item[PARENT_GENE] in gene_excluded:
                        gene_unclear[overlap_item[PARENT_GENE]] = dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]]
                        del gene_excluded[overlap_item[PARENT_GENE]]
                        # for tr in dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][TRANSCRIPT]:
                        #     del tr_excluded[tr]
                        # print("found single transcript: sitting on exluded gene " + overlap_item[PARENT_GENE] + " " + overlap_item[ID])
                        # gene_included[overlap_item[PARENT_GENE]] = dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]]
                        # del gene_excluded[overlap_item[PARENT_GENE]]
                        # for tr in dict_GE_TR[overlap_item[SOURCE]][overlap_item[PARENT_GENE]][TRANSCRIPT]: # remove all transcripts
                        #     if tr not in tr_included:
                        #         tr_included[tr] = dict_TR_GE[overlap_item[SOURCE]][overlap_item[ID]]
                        #     del tr_excluded[tr]
                    # else:
                        # print("found single transcript: sitting on included gene " + overlap_item[PARENT_GENE] + " " + overlap_item[ID])

    print("Genes included (before unclear): " + str(len(gene_included)))
    # try solve unclear
    remove_keys = []
    for gene in gene_unclear.keys():
        count_single = 0
        for tr in gene_unclear[gene][TRANSCRIPT]:
            if tr_excluded[tr] == "":
                count_single += 1
        if count_single >= 1: # the gene has at least 2 transcripts which are not overlapping
            gene_included[gene] = ''
            for tr in gene_unclear[gene][TRANSCRIPT]:
                tr_included[tr] = ''
                del tr_excluded[tr]
            remove_keys.append(gene)

    for k in remove_keys:
        del gene_unclear[k]

    fout = open("include_genes.txt","w")
    for g in gene_included:
        fout.write(g + "\n")
    fout.close()

    fout = open("exclude_genes.txt","w")
    for g in gene_excluded:
        fout.write(g + "\n")
    fout.close()

    fout = open("include_transcripts.txt","w")
    for g in tr_included:
        fout.write(g + "\n")
    fout.close()

    fout = open("exclude_transcripts.txt","w")
    for g in tr_excluded:
        fout.write(g + "\n")
    fout.close()

    with open("unclear.txt", 'w') as outfile:
        json.dump(gene_unclear, outfile, indent=4)

    # print("test")
    # for key in gene_included:
    #     if key in gene_excluded:
    #         print("doppelt include to exclude ge" + key)
    #
    # for key in tr_included:
    #     if key in tr_excluded:
    #         print("doppelt include to exclude tr" + key)
    #
    # for key in gene_excluded:
    #     if key in gene_included:
    #         print("doppelt exclude to include tr" + key)
    #
    # for key in tr_excluded:
    #     if key in tr_included:
    #         print("doppelt " + key)


    print("Genes included: " + str(len(gene_included)))
    print("Genes excluded: " + str(len(gene_excluded)))
    print("Transcripts included: " + str(len(tr_included)))
    print("Transcripts excluded: " + str(len(tr_excluded)))
    print(count_conflict)

    # # find genes which are both on included and excluded list
    # for key in gene_included:
    #     if key in gene_excluded:
    #         gene_unclear[key] = gene_excluded[key]
    #         gene_unclear[key][COUNT_EX] = 0
    #         gene_unclear[key][COUNT_IN] = 0
    #         for tr in gene_excluded[key][TRANSCRIPT]:
    #             if tr in tr_excluded:
    #                 gene_unclear[key][COUNT_EX] += 1
    #             else:
    #                 gene_unclear[key][COUNT_IN] += 1
    # print(len(gene_unclear))

if __name__ == "__main__":

    ################ Assgin unique keys #####################
    # asignUniqueKeys()
    #
    # insertKeys("mapping_keys.txt", "mapping_keys.json", "lncipedia_TR_keys.txt","TRANSCRIPT", "hg38", "TR")
    # insertKeys("mapping_keys.txt", "mapping_keys.json", "lncipedia_GE_keys.txt","GENE", "hg38", "GE")
    # insertKeys("mapping_keys.txt", "mapping_keys.json", "lncrnadb_TR_keys.txt","TRANSCRIPT", "hg38", "TR")
    # insertKeys("mapping_keys.txt", "mapping_keys.json", "lncrnadb_GE_keys.txt","GENE", "hg38", "GE")

    ################ Add gene entries #####################
    # addGeneEntrie("NONCODEv5_mouse_mm10_lncRNA.gtf")
    # addGeneEntrie("NONCODEv5_human_hg38_lncRNA.gtf")
    # addGeneEntrie("lncrnadb_mm10.gtf")
    # addGeneTranscriptEntry("lncipedia_5_0_hc_hg38_sorted.gtf")


   ################ Mapping keyes #####################

    # replaceKeyGTF("lncrnadb_mm10.gtf_withgenes",dict_mapping)
    # replaceKeyGTF("NONCODEv5_mouse_mm10_lncRNA.gtf_withgenes",dict_mapping)
    # replaceKeyGTF("NONCODEv5_human_hg38_lncRNA.gtf_withgenes",dict_mapping)

    # replaceKeyGTF("gencode.vM17.long_noncoding_RNAs.gtf",dict_mapping)
    # replaceKeyGTF("gencode.v28.long_noncoding_RNAs.gtf",dict_mapping)

    # replaceKeyGTF("hg38_long_noncoding_lncipedia.gtf",dict_mapping)
    # replaceKeyGTF("lncrnadb_hg38.gtf",dict_mapping)
    # replaceKeyGTF("lncrnadb_hg38.gtf","hg38")
    # replaceKeyGTF("hg38_noncode.gtf","hg38")
    # replaceKeyGTF("hg38_lncipedia.gtf","hg38")
    ################## Create big dict ###########################

    # build_nonredundantset("hg38")
    removeRedundancy("hg38")