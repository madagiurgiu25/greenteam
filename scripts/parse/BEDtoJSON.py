def parse(filein, fileout,classification,database,species):
    
    input=open(filein,"r")
    output=open(fileout,"w")
    header=["Chromosome","chromStart","chromEnd","Name","Score","Strand","thickStart","thickEnd","RGB","blockCount","blockSizes","blockStarts"]
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
        for i in range(0,len(entries)-1):
            output.write("\n\t\t\t\"Expression_"+header[i]+"\": \""+entries[i]+"\",")
        output.write("\n\t\t\t\"Expression_"+header[len(entries)-1].rstrip()+"\": \""+entries[len(entries)-1].rstrip()+"\"")
        output.write("\n\t\t}")  
    input.close()
    output.write("\n\t]")
    output.close()
    
dirpath="/home/offensperger/Schreibtisch/db/"
parse(dirpath+"NONCODEv5_hg38.lncAndGene.bed/data",dirpath+"pythontest.json","lncRNA","noncode","Homo sapiens")
