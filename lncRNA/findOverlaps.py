__author__ = 'Mada'

import os
import sys
from pathlib import Path
from pathlib import PurePosixPath
import json
import re

# one transcript_id
dict_lncRNA = {}

dict_visitedKeys = {}

# constants
GENCODE='gencode'
LNCIPEDIA='lncipedia'
NONCODE='noncode'
LNCRNADB='lncrnadb'
COUNT = 'count'
ALIAS='alias'
dict_index={GENCODE:2,LNCIPEDIA:2,NONCODE:2,LNCRNADB:2}
dict_index_gene={GENCODE:0,LNCIPEDIA:0,NONCODE:0,LNCRNADB:0}
dict_positions={GENCODE:0,LNCIPEDIA:1,NONCODE:2,LNCRNADB:3}
dict_positions_mouse={GENCODE:0,NONCODE:1,LNCRNADB:2}
SOURCE='source'


dict_db={ GENCODE: {},
    NONCODE: {},
    LNCIPEDIA: {},
    LNCRNADB: {}}

def readDictExons(filename):
    dict ={}
    with open(filename, "r") as f:
        for line in f:
            arr = line.strip().replace("\n","").split("\t")
            dict[arr[1]] = arr[0]
    return dict

def mapOverlaps(filename,primary_source):
    print("Start mapping " + primary_source)
    count = 0
    countall = 0
    with open(filename, "r") as f:
        for line in f:
            countall +=1
            
            arr = line.strip().replace("\n","").split("\t")
            countExons = arr[0]
            id1 = arr[1].strip()
            source = arr[2]
            id2 = arr[3].strip()

            _id1 = id1.split("@")[dict_index[primary_source]]
            _id2 = id2.split("@")[dict_index[source]]

            # print("Exons=" + countExons + "\t" + dict_db[primary_source][id1] + "\t" + _id1 + "\t" + dict_db[source][id2] + "\t" + _id2)

            # if transcript exons number match
            if id1 in dict_db[primary_source] and id2 in dict_db[source]:
                if countExons == dict_db[primary_source][id1] and countExons == dict_db[source][id2]:
                    # print("match")
                    count += 1
                    if _id1 not in dict_visitedKeys and _id2 not in dict_visitedKeys:
                        dict_lncRNA[_id1] = {}
                        dict_lncRNA[_id1][COUNT] = countExons
                        dict_visitedKeys[_id1] = ""
                        dict_lncRNA[_id1][ALIAS] = []
                        dict_lncRNA[_id1][SOURCE] = primary_source
                        dict_lncRNA[_id1][ALIAS].append((_id2,source))
                        dict_visitedKeys[_id2] = ""
                    elif _id1 not in dict_visitedKeys and _id2 in dict_visitedKeys:
                        if _id2 in dict_lncRNA:
                            dict_lncRNA[_id2][ALIAS].append((_id1,primary_source))
                            dict_visitedKeys[_id1] = ""
                        else:
                            # we have to find the key of _id2
                            for key in dict_lncRNA:
                                if (_id2,primary_source) in dict_lncRNA[key][ALIAS]:
                                    dict_lncRNA[key][ALIAS].append((_id1,primary_source))
                                    dict_visitedKeys[_id1] = ""
                                    break
                    elif _id1 in dict_visitedKeys and _id2 not in dict_visitedKeys:
                        if _id1 in dict_lncRNA:
                            dict_lncRNA[_id1][ALIAS].append((_id2,source))
                            dict_visitedKeys[_id2] = ""
                        else:
                            # we have to find the key of _id2
                            for key in dict_lncRNA:
                                if (_id1,primary_source) in dict_lncRNA[key][ALIAS]:
                                    dict_lncRNA[key][ALIAS].append((_id2,source))
                                    dict_visitedKeys[_id2] = ""
                                    break
            else:
                print("keys not found " + line)
    print("Overlaps that match also the number of transcripts are: " + str(count) + " out of " + str(countall) )

def convertDict2Matrix(fileout):

    fout = open(fileout,'w')
    fout2 =	open(fileout + "map",'w')
    fout.write("lncNamePrimary\texons_number\tGENCODE\tNONCODE\tLncipedia\tlncRNAdb\n")
    fout2.write("lncNamePrimary\texons_number\tGENCODE\tNONCODE\tLncipedia\tlncRNAdb\n")

    for key in dict_lncRNA:
        arr=[0,0,0,0]
        arr_string=["","","",""]
        arr[dict_positions[dict_lncRNA[key][SOURCE]]] = 1
        arr_string[dict_positions[dict_lncRNA[key][SOURCE]]] = key
        count_tr = dict_lncRNA[key][COUNT]
        for (id, source) in dict_lncRNA[key][ALIAS]:
            arr[dict_positions[source]] = 1
            arr_string[dict_positions[source]] = id
        fout.write(key + "\t" + str(count_tr) + "\t" + str(arr[0]) + "\t" + str(arr[1]) + "\t" + str(arr[2]) + "\t"+ str(arr[3]) + "\n")
        fout2.write(key + "\t" + str(count_tr) + "\t" + str(arr_string[0]) + "\t" + str(arr_string[1]) + "\t" + str(arr_string[2]) + "\t"+ str(arr_string[3]) +"\n")

    fout.close()
    fout2.close()

def convertDict2Matrix_mouse(fileout):

    fout = open(fileout,'w')
    fout2 =	open(fileout + "map",'w')
    fout.write("lncNamePrimary\texons_number\tGENCODE\tNONCODE\tlncRNAdb\n")
    fout2.write("lncNamePrimary\texons_number\tGENCODE\tNONCODE\tlncRNAdb\n")

    for key in dict_lncRNA:
        arr=[0,0,0]
        arr_string=["","",""]
        arr[dict_positions_mouse[dict_lncRNA[key][SOURCE]]] = 1
        arr_string[dict_positions_mouse[dict_lncRNA[key][SOURCE]]] = key
        count_tr = dict_lncRNA[key][COUNT]
        for (id, source) in dict_lncRNA[key][ALIAS]:
            arr[dict_positions_mouse[source]] = 1
            arr_string[dict_positions_mouse[source]] = id

        fout.write(key + "\t" + str(count_tr) + "\t" + str(arr[0]) + "\t" + str(arr[1]) + "\t" + str(arr[2]) + "\n")
        fout2.write(key + "\t" + str(count_tr) + "\t" + str(arr_string[0]) + "\t" + str(arr_string[1]) + "\t" + str(arr_string[2]) + "\n")
    fout.close()
    fout2.close()



def loadFiles_human():

	################# human

    # load the actual number of exons per transcript
    print("load dictionaries - numbers of exons per transcript per database source")
    dict_db[GENCODE] = readDictExons('gencode_hg38_short_exons.bed')
    dict_db[NONCODE] = readDictExons('noncode_hg38_short_exons.bed')
    dict_db[LNCRNADB] = readDictExons('lncrnadb_hg38_short_exons.bed')
    dict_db[LNCIPEDIA] = readDictExons('lncipedia_hg38_short_exons.bed')

    ##load overlaps for gencode
    print("map gencode and other sources")
    mapOverlaps('gencode_hg38_overlappAll_exons.bed',GENCODE)
    print("map noncode and other sources")
    mapOverlaps('noncode_hg38_overlappAll_exons.bed',NONCODE)
    print("map lncipedia and other sources")
    mapOverlaps('lncipedia_hg38_overlappAll_exons.bed',LNCIPEDIA)
    print("map lncrnadb and other sources")
    mapOverlaps('lncrnadb_hg38_overlappAll_exons.bed',LNCRNADB)

    # attach the rest of the
    print("attach the rest")
    for key in dict_db:
        for key_id in dict_db[key]:
            id  = key_id.split("@")[dict_index[key]]
            if id in dict_lncRNA:
                continue
            if id in dict_visitedKeys:
                continue

            dict_lncRNA[id] = {}
            dict_lncRNA[id][COUNT] = dict_db[key][key_id]
            dict_lncRNA[id][ALIAS] = []
            dict_lncRNA[id][SOURCE] = key

    print("convert to matrix and print to file")
    convertDict2Matrix("outOverlap.txt")

    print("convert to db json")
    dict_json=[]
    count=0
    for key in dict_lncRNA:
        dict_aux={}
        key_lnc = "LNC_HG38_" + str(count).zfill(6)
        count += 1
        dict_aux[key_lnc] = []
        dict_aux[key_lnc].append({SOURCE:dict_lncRNA[key][SOURCE],ALIAS:key})
        for (id,source) in dict_lncRNA[key][ALIAS]:
            dict_aux[key_lnc].append({SOURCE:source,ALIAS:id})
        dict_json.append(dict_aux)

    with open("overlap_hg38.json", 'w') as outfile:
        json.dump(dict_json, outfile, indent=4)
    outfile.close()

def loadFiles_mouse():

	# load the actual number of exons per transcript
    print("load dictionaries - numbers of exons per transcript per database source")
    dict_db[GENCODE] = readDictExons('gencode_mm10_short_exons.bed')
    dict_db[NONCODE] = readDictExons('noncode_mm10_short_exons.bed')
    dict_db[LNCRNADB] = readDictExons('lncrnadb_mm10_short_exons.bed')

    # load overlaps for gencode
    print("map gencode and other sources")
    mapOverlaps('gencode_mm10_overlappAll_exons.bed',GENCODE)
    print("map noncode and other sources")
    mapOverlaps('noncode_mm10_overlappAll_exons.bed',NONCODE)
    print("map lncrnadb and other sources")
    mapOverlaps('lncrnadb_mm10_overlappAll_exons.bed',LNCRNADB)

    # attach the rest of the
    print("attach the rest")
    for key in dict_db:
        for key_id in dict_db[key]:
            id  = key_id.split("@")[dict_index[key]]
            if id in dict_lncRNA:
                continue
            if id in dict_visitedKeys:
                continue

            dict_lncRNA[id] = {}
            dict_lncRNA[id][COUNT] = dict_db[key][key_id]
            dict_lncRNA[id][ALIAS] = []
            dict_lncRNA[id][SOURCE] = key

    print("convert to matrix and print to file")
    convertDict2Matrix_mouse("outOverlap_mm10.txt")

    print("convert to db json")
    dict_json=[]
    count=0
    for key in dict_lncRNA:
        dict_aux={}
        key_lnc = "LNC_MM10_" + str(count).zfill(6)
        count += 1
        dict_aux[key_lnc] = []
        dict_aux[key_lnc].append({SOURCE:dict_lncRNA[key][SOURCE],ALIAS:key})
        for (id,source) in dict_lncRNA[key][ALIAS]:
            dict_aux[key_lnc].append({SOURCE:source,ALIAS:id})
        dict_json.append(dict_aux)

    with open("overlap_mm10.json", 'w') as outfile:
        json.dump(dict_json, outfile, indent=4)
    outfile.close()


if __name__ == "__main__":
    loadFiles_human()
    # loadFiles_mouse()
