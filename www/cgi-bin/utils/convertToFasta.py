
def __getFasta__(fastaFile, keysFile,type):

    fin = open(keysFile,"r")
    dict = {}

    for l in fin.readlines():
        if l.startswith("#") == False:
            header = l.strip().split("\t")

            if type == 'transcript':
                dict[header[3]] = ''

    fin.close()
    print(dict)

    fin = open(fastaFile,'r')
    fout = open("out.fasta",'w')

    continueRead = False
    for l in fin.readlines():

        if l.startswith(">"):
            if continueRead == True:
                continueRead = False
            header = l[1:].split("|")
            if header[0] in dict:
                fout.write(l)
                continueRead = True
        else:
            if continueRead == True:
                fout.write(l)

    fin.close()
    fout.close()

if __name__ == "__main__":
    # __getFasta__("gencode.vM17.transcripts.fa","mm10.bed","transcript")
    __getFasta__("gencode.v28.transcripts.fa", "hg38.bed", "transcript")