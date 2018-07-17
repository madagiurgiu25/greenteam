__author__ = 'Mada'
# !/usr/bin/python3
# !C:\\Python34\\python.exe

import computeCoverage
import sys

if __name__ == "__main__":
    if len(sys.argv) == 3:
        # compute(json_species,gene_abundance_sample)
        computeCoverage.computePerSample(sys.argv[2],sys.argv[1])