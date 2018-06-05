__author__ = 'Mada'
#!/usr/bin/python3
#!C:\\Python34\\python.exe

########## find overlaps between genes
import json
import intervaltree
import re
import operator

# headers in GTF file
gtfdict = {'gene_id': "", 'transcript_id': "", 'transcript_type': "", 'transcript_name': "", 'gene_name': "",
           'gene_type': "", 'transcript_type':"", 'tag':"", 'exon_number':"", 'Parent':"", 'level':""}

dict_words = {'geneID': "gene_id", 'transcriptID': "transcript_id",'exonnumber':"exon_number"}

GENE_ID = 'geneID'
TRANSCRIPT_ID = 'transcriptID'
EXON_NR = 'exonnumber'

EXON = 'exon'
GENE = 'gene'
TRANSCRIPT = 'transcript'
UTR = 'UTR'
CDS = 'CDS'

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# require json files
def findOverlaps(file1, file2=None):

    count = 0
    count_overlaps = 0
    visited=[]
    fout=open("overlaps", "w")
    # fout.write("gene_id_1\tchr_1\tstart_1\tstop_1\tstrand_1\tgene_id_2\tchr_2\tstart_2\tstop_2\tstrand_2\tcount_overlap\tpercent_1\tpercent_2")
    fout.write("gene_id_1\tstrand_1\tgene_id_2\tstrand_2\tcount_overlap\tfraction_1\tfraction_2")
    if file2 == None: # search overlaps between genes in the same file
        try:
            with open(file1, 'r') as f:
                dict_genes = json.load(f)
                print("Total number of genes: " +  str(len(dict_genes)))
                for index, gene in enumerate(dict_genes):
                    print(count)
                    count = count + 1
                    gene_1=gene["gene_id"]
                    visited.append(gene_1)
                    for index_2, gene_aux in enumerate(dict_genes[index+1:]):
                        gene_2=gene_aux["gene_id"]
                        if gene["chr"] == gene_aux["chr"]:
                            overlap = getOverlap([int(gene["start"]),int(gene["stop"])],[int(gene_aux["start"]),int(gene_aux["stop"])])
                            if overlap > 0:
                                print(gene_1 + "\t" + gene_2 + "\t" + str(overlap))
                                gene_1_overlap=round((float)(overlap+1)/(int(gene["stop"])-int(gene["start"]) +1),4)
                                gene_2_overlap=round((float)(overlap+1)/(int(gene_aux["stop"])-int(gene_aux["start"]) +1),4)
                                overlap_len = overlap +1
                                fout.write(gene_1 + "\t" + gene['strand'] + "\t" + gene_2 + "\t" + gene_aux["strand"] + "\t" + str(overlap_len) + "\t" + str(gene_1_overlap) + "\t" + str(gene_2_overlap) + "\n")
                                count_overlaps += 1
                                # fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n").format(gene_1,
                                #                                                gene["strand"],
                                #                                                gene_2,
                                #                                                gene_aux["strand"],
                                #                                                str(overlap),
                                #                                                str(gene_1_overlap),
                                #                                                str(gene_2_overlap)
                                #                                                )

                print("Total number of genes: " +  str(len(dict_genes)))
                print("Number of overlaps: " + str(count_overlaps))
        except FileNotFoundError:
            print("your file does not exist")

def findOverlaps_lnc2gene(file):

    # construct interval tree
    it_regions  = intervaltree.IntervalTree()
    try:
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                arr = line.strip().replace("\n","").split("\t")
                chr=arr[0]
                start=int(arr[3])
                stop=int(arr[4])
                type=arr[2]

                geneID=re.search("gene_id \"(.*?)\";", arr[8], flags=0)
                transcriptID=re.search("transcript_id \"(.*?)\";", arr[8], flags=0)
                exonnumber=re.search("exon_number (\d);",arr[8],flags=0)

                if stop > start:
                    dict = {}
                    dict['chr'] = chr
                    dict['type'] = type
                    if geneID is not None:
                        dict['geneID'] = geneID.group(1)
                    if transcriptID is not None:
                        dict['transcriptID'] = transcriptID.group(1)
                    if exonnumber is not None:
                        dict['exonnumber'] = exonnumber.group(1)

                    it_regions[start:stop] = dict
                    # print(it_regions[start:stop])
                   # print(it_regions[start:stop])

        dict={}
        ##### find overlaps
        for (start,stop,data) in it_regions:
            if data['type'] == 'gene' and GENE_ID in dict_words and len(data[GENE_ID]) > 3 and data[GENE_ID][0:3] == "LNC":
                fi = it_regions.search(start,stop,{'chr':data['chr']})
                arr = {'exon':0,'CDS':0,'UTR':0}
                for (rstart,rstop,rdata) in fi:
                    if rdata[GENE_ID] != data[GENE_ID] and rdata[GENE_ID][0:3] == "LNC":
                        if rdata['type'] == EXON:
                            arr[EXON] += (rstop-rstart)/(stop-start)
                        elif rdata['type'] == UTR:
                            arr[UTR] += (rstop - rstart) / (stop - start)
                        elif rdata['type'] == CDS:
                            arr[CDS] += (rstop - rstart) / (stop - start)
                dict[data[GENE_ID]] = {}
                dict[data[GENE_ID]]['arr'] = arr

        fout = open("test_chr18.txt", 'w')
        # find max
        for item in dict:
            dict[item]['HIT'] = max(dict[item]['arr'].items(), key=operator.itemgetter(1))[0]
            dict[item]['SCORE'] = dict[item]['arr'][dict[item]['HIT']]
            if dict[item]['SCORE'] == 0:
                fout.write(item + "\t" + "INTERGENIC" + "\t" + str(dict[item]['SCORE']) + "\n")
            else:
                fout.write(item + "\t" + dict[item]['HIT']  + "\t" + str(dict[item]['SCORE']) + "\n")

        with open("test_chr18.json", 'w') as outfile:
            json.dump(dict, outfile, indent=4)
        outfile.close()

    except FileNotFoundError:
        print("your file does not exist")

if __name__ == "__main__":
    findOverlaps("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/mapping/mm10_primary_assembly_and_lncRNA.gtf")
    # findOverlaps("mm10_primary_assembly_and_lncRNA.json")
    # findOverlaps_lnc2gene("chr18.gtf")