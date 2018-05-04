#!/usr/bin/python3
#!C:\\Python34\\python.exe

import os
import sys
from pathlib import Path
from pathlib import PurePosixPath
import json

# path to the current script
path_to_script = os.path.realpath(__file__)

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

    except FileNotFoundError:
        print("your file does not exist")

if __name__ == "__main__":
    file = sys.argv[1]
    parseBed2Json(file)
