#!/usr/bin/python3
# !C:\\Python34\\python.exe

import os
import sys
from pathlib import Path
from pathlib import PurePosixPath
import json
import re

# path to the current script
path_to_script = os.path.realpath(__file__)

# Entry type in gtf
EXON = 'exon'
GENE = 'gene'
TRANSCRIPT = 'transcript'

# dictionary constants
START = 'start'
STOP = 'stop'
CHR = 'chr'
STRAND = 'strand'
GENE_ID = 'gene_id'
TRANSCRIPT_ID = 'transcript_id'
TRANSCRIPTS = 'transcripts_list'
EXONS = 'exon_list'
TAG = 'tag'
GENE_TYPE = 'gene_type'
GENE_NAME = 'gene_name'
TRANSCRIPT_NAME = 'transcript_name'
TRANSCRIPT_TYPE = 'transcript_type'

# headers in GTF file
gtfdict = {'gene_id': "", 'transcript_id': "", 'transcript_type': "", 'transcript_name': "", 'gene_name': "",
           'gene_type': "", 'transcript_type':"", 'tag':"", 'exon_number':""}

# headers in GTF file for GENE entry
gtfdict_gene = {'gene_id': "", 'gene_name': "", 'gene_type': "", 'tag':""}

# headers in GTF file for TRANSCRIPT entry
gtfdict_transcript = {'transcript_id': "", 'transcript_type': "", 'transcript_name': "", 'transcript_type':"", 'tag':""}

# headers in GTF file for EXON entry
gtfdict_exon = {'exon_number': "", 'tag':""}


def infoString2Dict(infoString):

    # meta info
    info = infoString.split(";") # Ex: gene_id "ENSG00000243485.5_4"; gene_type "lincRNA"; gene_name "RP11-34P13.3"; level 2; tag "ncRNA_host";
    dictInfo={}
    for i in info:

        items = re.split(" |=",i.strip())

        # only if we have 'key = value' string / overcome empty strings
        if len(items) == 2:
            (key,value) = items
            if key in gtfdict:
                dictInfo[key] = value.strip().replace("\"","")

    return dictInfo

def parseRow2Exon(arr):

    dict = {}
    dict[CHR] = arr[0]
    dict[START] = arr[3]
    dict[STOP] = arr[4]
    dict[STRAND] = arr[6]

    # meta info
    info = arr[8].split(";") # Ex: gene_id "ENSG00000243485.5_4"; gene_type "lincRNA"; gene_name "RP11-34P13.3"; level 2; tag "ncRNA_host";

    for i in info:

        items = re.split(" |=",i.strip())
        # only if we have 'key = value' string / overcome empty strings
        if len(items) == 2:
            (key,value) = items
            if key in gtfdict_exon:
                dict[key] = value.strip().replace("\"","")

    return dict

def parseRow2Transcript(arr):

    dict = {}
    dict[CHR] = arr[0]
    dict[START] = arr[3]
    dict[STOP] = arr[4]
    dict[STRAND] = arr[6]
    dict[EXONS] = []

    # meta info
    info = arr[8].split(";") # Ex: gene_id "ENSG00000243485.5_4"; gene_type "lincRNA"; gene_name "RP11-34P13.3"; level 2; tag "ncRNA_host";

    for i in info:

        items = re.split(" |=",i.strip())
        # only if we have 'key = value' string / overcome empty strings
        if len(items) == 2:
            (key,value) = items
            if key in gtfdict_transcript:
                dict[key] = value.strip().replace("\"","")

    return dict

def parseRow2Gene(arr):

    dict = {}
    dict[CHR] = arr[0]
    dict[START] = arr[3]
    dict[STOP] = arr[4]
    dict[STRAND] = arr[6]
    dict[TRANSCRIPTS] = []

    # meta info
    info = arr[8].split(";") # Ex: gene_id "ENSG00000243485.5_4"; gene_type "lincRNA"; gene_name "RP11-34P13.3"; level 2; tag "ncRNA_host";

    for i in info:

        items = re.split(" |=",i.strip())

        # only if we have 'key = value' string / overcome empty strings
        if len(items) == 2:
            (key,value) = items
            if key in gtfdict_gene:
                dict[key] = value.strip().replace("\"","")

    return dict

# extract from gft file format entry of gene_type x (miRNA, protein_coding, lincRNA)
def parseGTF(gtffile,bedoutput,jsonoutput, assembly, source, type=None):

    bigDict = []

    genesAlreadyParsed = {}
    transcriptsAlreadyParsed = {}

    # test if the file exitst
    try:
        with open(gtffile, "r") as f:
            for line in f:
                if line.startswith("#") is False: # ignore header lines
                    arr = line.strip("\n").split("\t")

                    # parse the info string in the entry
                    infoDict = infoString2Dict(arr[8])

                    # filter for gene type
                    if type is not None and GENE_TYPE not in infoDict:
                        continue

                    if type is not None and infoDict[GENE_TYPE] != type:
                        continue

                    ############### Update information for the parents already inserted in the bigDict
                    # GENE
                    if arr[2] == GENE and infoDict[GENE_ID] in genesAlreadyParsed:
                        # insert chr / start / stop / strand
                        for item in bigDict:
                            if item[GENE_ID] == infoDict[GENE_ID]:
                                item[CHR] = arr[0]
                                item[START] = arr[3]
                                item[STOP] = arr[4]
                                item[STRAND] = arr[6]
                                break
                        continue
                    # TRANSCRIPT
                    if arr[2] == TRANSCRIPT and infoDict[TRANSCRIPT_ID] in transcriptsAlreadyParsed:
                        # insert chr / start / stop / strand
                        for item in bigDict:
                            if item[GENE_ID] == infoDict[GENE_ID]:
                                for t in item[TRANSCRIPTS]:
                                    if t[TRANSCRIPT_ID] == infoDict[TRANSCRIPT_ID]:
                                        t[CHR] = arr[0]
                                        t[START] = arr[3]
                                        t[STOP] = arr[4]
                                        t[STRAND] = arr[6]
                                        break
                        continue

                    ############## Add element (Gene/Transcript/Exon
                    # GENE
                    if arr[2] == GENE: # gene can contain multiple transcripts
                        newGene = parseRow2Gene(arr)
                        bigDict.append(newGene)
                        genesAlreadyParsed[newGene[GENE_ID]] = ''
                        continue

                    # TRANSCRIPT
                    if arr[2] == TRANSCRIPT:
                        newTranscript = parseRow2Transcript(arr)

                        # Parent GENE already in bigDict
                        if infoDict[GENE_ID] in genesAlreadyParsed:
                            for item in bigDict:
                                if item[GENE_ID] == infoDict[GENE_ID]:
                                    item[TRANSCRIPTS].append(newTranscript)
                                    transcriptsAlreadyParsed[newTranscript[TRANSCRIPT_ID]] = ''
                                    break
                        else:
                            # create and parent gene
                            parentGene = {}
                            parentGene[GENE_ID] = infoDict[GENE_ID]
                            parentGene[GENE_TYPE] = infoDict[GENE_TYPE]
                            parentGene[GENE_NAME] = infoDict[GENE_NAME]
                            parentGene[TRANSCRIPTS] = []

                            bigDict.append(parentGene)
                            genesAlreadyParsed[parentGene[GENE_ID]] = ''

                            # add transcript
                            parentGene[TRANSCRIPTS].append(newTranscript)
                            transcriptsAlreadyParsed[newTranscript[TRANSCRIPT_ID]] = ''
                        continue

                    # EXON
                    if arr[2] == EXON:

                        newExon = parseRow2Exon(arr)
                        # Check if parents (GENE and TRANSCRIPT) are already inserted
                        if infoDict[GENE_ID] in genesAlreadyParsed and infoDict[TRANSCRIPT_ID] in transcriptsAlreadyParsed:
                            # Search element
                            for item in bigDict:
                                if item[GENE_ID] == infoDict[GENE_ID]:
                                    for t in item[TRANSCRIPTS]:
                                        if t[TRANSCRIPT_ID] == infoDict[TRANSCRIPT_ID]:
                                            t[EXONS].append(newExon)
                                            break

                        elif infoDict[GENE_ID] in genesAlreadyParsed and infoDict[TRANSCRIPT_ID] not in transcriptsAlreadyParsed:
                            # GENE parent exists but not the Transcript
                            # create and parent transcript
                            parentTranscript = {}
                            parentTranscript[TRANSCRIPT_ID] = parentTranscript[TRANSCRIPT_ID]
                            parentTranscript[TRANSCRIPT_NAME] = parentTranscript[TRANSCRIPT_NAME]
                            parentTranscript[EXONS] = []

                            transcriptsAlreadyParsed[parentTranscript[TRANSCRIPT_ID]] = ''

                            # add the exon to transcript
                            parentTranscript[EXONS].append(newExon)

                            # add the transcript to parent gene
                            for item in bigDict:
                                if item[GENE_ID] == infoDict[GENE_ID]:
                                    item[TRANSCRIPTS].append(parentTranscript)

                        elif infoDict[GENE_ID] not in genesAlreadyParsed:
                            # nor GENE nor TRANSCRIPT PArents are available

                            # create parent transcript
                            parentTranscript = {}
                            parentTranscript[TRANSCRIPT_ID] = parentTranscript[TRANSCRIPT_ID]
                            parentTranscript[TRANSCRIPT_NAME] = parentTranscript[TRANSCRIPT_NAME]
                            parentTranscript[EXONS] = []
                            transcriptsAlreadyParsed[parentTranscript[TRANSCRIPT_ID]] = ''

                            # add the exon to transcript
                            parentTranscript[EXONS].append(newExon)

                            # create and parent gene for transcript
                            parentGene = {}
                            parentGene[GENE_ID] = infoDict[GENE_ID]
                            parentGene[GENE_TYPE] = infoDict[GENE_TYPE]
                            parentGene[GENE_NAME] = infoDict[GENE_NAME]
                            parentGene[TRANSCRIPTS] = []

                            # add transcript to parent gene
                            parentGene[TRANSCRIPTS].append(parentTranscript)

                            bigDict.append(parentGene)
                            genesAlreadyParsed[parentGene[GENE_ID]] = ''
                        continue
        f.close()

        ############ Parse and write data to bed
        parseGTF2BED(bigDict,bedoutput, assembly, source)

        ############ DumpJson
        with open(jsonoutput, 'w') as outfile:
            json.dump(bigDict, outfile, indent=4)
        outfile.close()

    except FileNotFoundError:
        print("your file does not exist")


def parseGTF2BED(gtfjson,bedoutput, assembly, source):

    try:

        fout = open(bedoutput,'w')
        # write header
        fout.write("#chr\tstart_exon\tstop_exon\tstrand\tsource\tassembly\tstart_gene\tstop_gene\tgene_id\tgene_name\tgene_type\tstart_transcript\tstop_transcript\ttranscript_name\ttranscript_type\n")

        for item in gtfjson:
            for t in item[TRANSCRIPTS]:
                for exon in t[EXONS]:
                    toWrite = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n"\
                        .format(exon[CHR],
                                exon[START],
                                exon[STOP],
                                exon[STRAND],
                                source,
                                assembly,
                                item[START],
                                item[STOP],
                                item[GENE_ID],
                                item[GENE_NAME],
                                item[GENE_TYPE],
                                t[START],
                                t[STOP],
                                t[TRANSCRIPT_ID],
                                t[TRANSCRIPT_NAME],
                                t[TRANSCRIPT_TYPE])
                    fout.write(toWrite)
        fout.close()
    except FileNotFoundError:
        print("your file does not exist")


def parseBed2Json(bedfile):
    # test if the file exitst
    try:
        # try resolve path and get the absolute path
        my_abs_path_file = Path(bedfile).resolve()

        ###### Open bed file and parse to json #####
        jsonData = []

        # open bed file
        infile = open(bedfile, 'r').readlines()

        # read header
        header = infile.pop(0).strip('\n').split("\t")
        for line in infile:

            # get row and convert to array
            data = line.strip('\n').split("\t")

            # test if the row is complete and have the same size as the header
            if len(data) == len(header):
                datadict = {}

                # iterate over my bedfile row
                for i, element in enumerate(data):
                    datadict[header[i]] = element

                # push into array the dictionary structure
                jsonData.append(datadict)

        ######## write into file the json structure

        # bed filename (without path)
        filename = PurePosixPath(my_abs_path_file).name

        # json file name
        jsonFile = filename + ".json"
        # test if file ends with .bed
        if filename.endswith('.bed'):
            jsonFile = filename.split(".bed")[0] + ".json"

        with open(jsonFile, 'w') as outfile:
            json.dump(jsonData, outfile, indent=4)
        outfile.close()

    except FileNotFoundError:
        print("your file does not exist")


if __name__ == "__main__":
    file = sys.argv[1]

    ############# parse generic BED file with header to json
    # parseBed2Json(file)

    file = sys.argv[1]
    assembly = sys.argv[2]
    source = sys.argv[3]
    bedfile = sys.argv[4]
    jsonfile = sys.argv[5]
    type = sys.argv[6]

    ########### parse GTF file to BED and json
    parseGTF(file,bedfile,jsonfile,assembly,source,type)
