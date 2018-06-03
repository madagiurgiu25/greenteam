__author__ = 'Mada'
#!/usr/bin/python3
#!C:\\Python34\\python.exe

########## find overlaps between genes
import json

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# require json files
def findOverlaps(file1, file2=None):

    count = 0
    visited=[]
    fout=open("overlaps", "w")
    # fout.write("gene_id_1\tchr_1\tstart_1\tstop_1\tstrand_1\tgene_id_2\tchr_2\tstart_2\tstop_2\tstrand_2\tcount_overlap\tpercent_1\tpercent_2")
    fout.write("gene_id_1\tstrand_1\tgene_id_2\tstrand_2\tcount_overlap\tfraction_1\tfraction_2")
    if file2 == None: # search overlaps between genes in the same file
        try:
            with open(file1, 'r') as f:
                dict_genes = json.load(f)
                for gene_1 in dict_genes:
                    print(count)
                    count = count + 1
                    visited.append(gene_1)
                    for gene_2 in dict_genes:
                        if gene_1 != gene_2 and gene_2 not in visited:
                            if dict_genes[gene_1]["chr"] == dict_genes[gene_2]["chr"]:
                                overlap = getOverlap([int(dict_genes[gene_1]["start"]),int(dict_genes[gene_1]["stop"])],[int(dict_genes[gene_2]["start"]),int(dict_genes[gene_2]["stop"])]
                                if overlap > 0:
                                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n").format(gene_1,
                                                                                   dict_genes[gene_1]["strand"],
                                                                                   gene_2,
                                                                                   dict_genes[gene_1]["strand"],
                                                                                   str(overlap+1),
                                                                                   str((float)(overlap+1)/(dict_genes[gene_1]["stop"]-dict_genes[gene_1]["start"] +1)),
                                                                                   str((float)(overlap+1)/(dict_genes[gene_2]["stop"]-dict_genes[gene_2]["start"] +1))
                                                                                   )


        except FileNotFoundError:
            print("your file does not exist")

if __name__ == "__main__":
    findOverlaps("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/data/mapping/mm10_primary_assembly_and_lncRNA.json")