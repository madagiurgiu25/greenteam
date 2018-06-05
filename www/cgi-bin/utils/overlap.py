__author__ = 'Mada'
#!/usr/bin/python3
#!C:\\Python34\\python.exe

########## find overlaps between genes
import json
import intervaltree
import re

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
                geneID=re.match("gene_id \"(\w+)\";", arr[8], flags=0)
                transcriptID=re.match("transcript_id \"(\w+)\";", arr[8], flags=0)
                exonnumber=re.match("exon number (\d);",arr[8],flags=0)

                if stop > start:
                    dict = {}
                    dict['chr'] = chr
                    dict['type'] = type
                    if geneID is not None:
                        dict['geneUD'] = geneID
                    if transcriptID is not None:
                        dict['transcriptID'] = transcriptID
                    if exonnumber is not None:
                        dict['exonnumber'] = exonnumber

                    it_regions[start:stop] = dict
                    # print(it_regions[start:stop])
                   # print(it_regions[start:stop])


        ##### find overlaps

        for (start,stop,data) in it_regions:
            # print(p)
            if data["geneID"].startWith("LNC"):
                print("do")
                # fi = find("LNC")
        # tree.intersect( start, end, lambda x: out.append(x) )
        # list = [ (x.start, x.end) for x in out ]


    except FileNotFoundError:
        print("your file does not exist")

if __name__ == "__main__":
    # findOverlaps("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/mapping/mm10_primary_assembly_and_lncRNA.json")
    # findOverlaps("mm10_primary_assembly_and_lncRNA.json")
    findOverlaps_lnc2gene("chr18.gtf")