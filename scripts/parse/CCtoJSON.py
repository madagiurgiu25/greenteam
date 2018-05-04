
def parse(filein, fileout,classification,database,species):
    
    input=open(filein,"r")
    output=open(fileout,"w")
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
        cclass=''
        if (entries[1][0]=='1') :
            cclass+="Antisense"
        if(entries[1][1]=='1'):
            if(cclass!=''): 
                cclass+=","
            cclass+="Exonic"
        if(entries[1][2]=='1'):
            if(cclass!=''):
                cclass+=","
            cclass+="Sense_No_Exonic"
        if(entries[1][3]=='1'):
            if(cclass!=''):
                cclass+=","
            cclass+="Intergenic"
        output.write("\n\t\t\t\"Class\": \""+cclass+"\"\n\t\t}")
    input.close()
    output.write("\n\t]")
    output.close()
    
dirpath="/home/offensperger/Schreibtisch/db/"
parse(dirpath+"NONCODEv5_mouse_gene_cc/data",dirpath+"NONCODEv5_mouse_gene_cc.json","lncRNA","noncode","Mus musculus")
parse(dirpath+"NONCODEv5_human_gene_cc/data",dirpath+"NONCODEv5_human_gene_cc.json","lncRNA","noncode","Homo sapiens")