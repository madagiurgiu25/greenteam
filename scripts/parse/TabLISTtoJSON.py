
def parse(filein, fileout,classification,database,species,annotation):
    
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
        if (database=="circ2Traits"):
             output.write("\n\t\t\t\""+annotation+"\": \""+str(entries[1].rstrip().substring(0, entries[1].length()-1).split(","))+"\"\n\t\t}")
        else:
             output.write("\n\t\t\t\""+annotation+"\": \""+str(entries[1].rstrip().split(","))+"\"\n\t\t}")
    input.close()
    output.write("\n\t]")
    output.close()
    
dirpath="/home/offensperger/Schreibtisch/db/"
parse(dirpath+"NONCODEv5_human.func/data",dirpath+"NONCODEv5_human.func.json","lncRNA","noncode","Homo sapiens","Function")
parse(dirpath+"NONCODEv5_mouse.func/data",dirpath+"NONCODEv5_mouse.func.json","lncRNA","noncode","Mus musculus","Function")
parse(dirpath+"all_disease_set.txt",dirpath+"all_disease_set.json","miRNA","circ2Traits","Homo sapiens","Disease_targets")