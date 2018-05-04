def parse(filein, fileout,classification,header,database,species):
    
    input=open(filein,"r")
    output=open(fileout,"w")
    output.write("[")
    started=False
    addtoseq=True
    sequence=''
    for line in input:
       
        if line.startswith(">"):
            if sequence !='':
                output.write("\n\t\t\t\"Sequence\": \""+sequence.upper()+"\"\n\t\t}")
                sequence=''
            if header==0: 
                if started:
                    output.write(",")
                else:
                    started=True
                output.write("\n\t\t{")
                output.write("\n\t\t\t\"Type\": \""+classification+"\",")
                output.write("\n\t\t\t\"Database\": \""+database+"\",")
                output.write("\n\t\t\t\"Species\": \""+species+"\",")
                output.write("\n\t\t\t\"Name\": \""+line[1:].rstrip()+"\",")
            if header==1:
                split=line[1:].split()
                if((split[2]=="Homo" and split[3]=="sapiens") or (split[2]=="Mus" and split[3]=="musculus")):
                    addtoseq=True
                    if(started):
                        output.write(",")
                    started=True
                    output.write("\n\t\t{")
                    output.write("\n\t\t\t\"Database\": \""+database+"\",")
                    output.write("\n\t\t\t\"miRBase_id\": \""+split[1]+"\",")
                    output.write("\n\t\t\t\"Species\": \""+split[2]+" "+split[3]+"\",")
                    output.write("\n\t\t\t\"Name\": \""+split[0]+"\",")
                else:
                    addtoseq=False
        else:
            if addtoseq:
                sequence+=line.rstrip()
    input.close()
    if sequence!='':
        output.write("\n\t\t\t\"Sequence\": \""+sequence.upper()+"\"\n\t\t}")
        sequence=""
    output.write("\n]")
    output.close()
    
dirpath="/home/offensperger/Schreibtisch/db/"
parse(dirpath+"mature.fa",dirpath+"pythontest.json","lncRNA",1,"lncipedia","Homo sapiens")
parse(dirpath+"lncipedia_5_0.fasta",dirpath+"lncipedia_5_0.json","lncRNA",0,"lncipedia","Homo sapiens")
parse(dirpath+"mature.fa",dirpath+"mature.json","miRNA",1,"miRBase",'')
parse(dirpath+"NONCODEv5_mouse.fa",dirpath+"NONCODEv5_mouse.json","lncRNA",0,"noncode","Mus musculus")
parse(dirpath+"NONCODEv5_human.fa",dirpath+"NONCODEv5_human.json","lncRNA",0,"noncode","Homo sapiens")
