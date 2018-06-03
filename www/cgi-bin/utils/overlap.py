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
                for gene in dict_genes:
                    print(count)
                    count = count + 1
                    gene_1=gene["gene_id"]
                    visited.append(gene_1)
                    for gene_aux in dict_genes:
                        gene_2=gene_aux["gene_id"]
                        if gene_1 != gene_2 and gene_2 not in visited:
                            if gene["chr"] == gene_aux["chr"]:
                                overlap = getOverlap([int(gene["start"]),int(gene["stop"])],[int(gene_aux["start"]),int(gene_aux["stop"])])
                                if overlap > 0:
                                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n").format(gene_1,
                                                                                   gene["strand"],
                                                                                   gene_2,
                                                                                   gene_aux["strand"],
                                                                                   str(overlap+1),
                                                                                   str((float)(overlap+1)/(gene["stop"]-gene["start"] +1)),
                                                                                   str((float)(overlap+1)/(gene_aux["stop"]-gene_aux["start"] +1))
                                                                                   )


        except FileNotFoundError:
            print("your file does not exist")

if __name__ == "__main__":
    findOverlaps("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/mapping/mm10_primary_assembly_and_lncRNA.json")