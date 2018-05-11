
def parse(filein, fileout,classification,database,species):
    
    input=open(filein,"r")
    output=open(fileout,"w")
    print(input.readline)
    header=input.readline().split("\t")
    output.write("[")
    started=False
    for line in input:
        if started:
            output.write(",")
        else: 
            started=True 
        output.write("\n\t\t{")
        output.write("\n\t\t\t\"Database\": \""+database+"\",")
        output.write("\n\t\t\t\"Species\": \""+species+"\",")
        output.write("\n\t\t\t\"Type\": \""+classification+"\",")
        entries=line.split("\t")
        output.write("\n\t\t\t\"Name\": \""+entries[0]+"\",")
        for i in range(1,len(entries)-1):
            output.write("\n\t\t\t\"Expression_"+header[i]+"\": \""+entries[i]+"\",")
        output.write("\n\t\t\t\"Expression_"+header[len(entries)-1].rstrip()+"\": \""+entries[len(entries)-1].rstrip()+"\"")
        output.write("\n\t\t}")  
    input.close()
    output.write("\n\t]")
    output.close()
    
dirpath="/home/offensperger/Schreibtisch/db/"

parse(dirpath+"E-MTAB-1969.processed.1/normalized_data.txt",dirpath+"E-MTAB-1969.processed.1/normalized_data.json","miRNA","ArrayExpress","Homo sapiens")

parse(dirpath+"NONCODEv5_human_exosome_gene.fpkm/data",dirpath+"NONCODEv5_human_exosome_gene.fpkm.json","Gene","noncode","Homo sapiens")
parse(dirpath+"NONCODEv5_human_exosome_lncRNA.fpkm/data",dirpath+"NONCODEv5_human_exosome_lncRNA.fpkm.json","lncRNA","noncode","Homo sapiens")

parse(dirpath+"NONCODEv5_mouse.gene.exp/data",dirpath+"NONCODEv5_mouse.gene.exp.json","Gene","noncode","Mus musculus")
parse(dirpath+"NONCODEv5_mouse.lncRNA.exp/data",dirpath+"NONCODEv5_mouse.lncRNA.exp.json","lncRNA","noncode","Mus musculus")

parse(dirpath+"NONCODEv5_human.gene.exp/data",dirpath+"NONCODEv5_human.gene.exp.json","Gene","noncode","Homo sapiens")
parse(dirpath+"NONCODEv5_human.lncRNA.exp/data",dirpath+"NONCODEv5_human.lncRNA.exp.json","lncRNA","noncode","Homo sapiens")