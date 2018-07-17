#!/usr/bin/env Rscript
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
library(plyr)
library(ggfortify)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyr)
library(dplyr)
library("NOISeq")
library(stringr)


setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/greenteam/scripts/diffexp/")
source("filtering.R", chdir = T)
source("plotsNormalization.R", chdir = T)
source("dataTransformation_library.R", chdir = T)


args = commandArgs(trailingOnly = TRUE)
if (length(args) <= 0) {
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
}

sample = args[1]
subset_size = as.numeric(args[2])
data_directory = (args[3])
working_dir = (args[4])
sampleNames_conditions = data.frame(read.csv2(file = args[5],sep = "\t",header = F))
anno =  data.frame(read.csv2(file = args[6], sep="\t", header = F))


print(sample)
print(subset_size)
print(data_directory)
print(working_dir)
print("Sample_Conditions:" + args[5])
print("Anno file:" + args[6])



#sample = "Mlet7"
#subset_size = 5
#data_directory = (
#  '/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/StringTie'
#)
#working_dir = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/DiffExp/Mlet7/NOISeq')
#sampleNames_conditions = data.frame(
#  read.csv2(file = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/samples/mlet7.txt", sep = "\t", header = F)
#)
#anno= data.frame(read.csv2(file = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/anno/mm10_primary_assembly_and_lncRNA_supergenes_exclude20.gtf_noiseq", sep="\t", header = F))
#size = 5

setwd(working_dir)

print("read input...")
sampleNames_conditions$V5 <- as.numeric(as.character(sampleNames_conditions$V5))
df_config <- sampleNames_conditions[order(sampleNames_conditions$V5), ]
df_config <- df_config[,3:6]
sampleConditions <- unique(data.frame(V5=as.numeric(as.character(df_config$V5)), V4=as.factor(df_config$V4)))
conditions <- dim(sampleConditions)[1]
size <- as.numeric(dim(unique(data.frame(df_config$V6)))[1])

list_sampleNames <- data.frame(aggregate(V3 ~ V4, data = df_config, paste, collapse = ","))
print(list_sampleNames)
list_replicates <- data.frame(aggregate(V6 ~ V4, data = df_config, paste, collapse = ","))
print(list_replicates)
list_sampleNames <- merge(list_sampleNames,sampleConditions)
list_sampleNames <- merge(list_sampleNames,list_replicates)

print(list_sampleNames)

sampleList <- combn(size, subset_size, simplify = FALSE)

for (i in 1:length(sampleList)) {
  
  # get all files
  vector <- unlist(sampleList[i], use.names = F)
  print("process subset....")
  print(vector)
  
  samplePaths <- vector()
  samplePaths_TPM <- vector()
  groupconditions <- vector()
  replicates <- vector()
  
  for (j in 1:dim(list_sampleNames)[1]){
    
    ########################################## Format input #####################################################
    
    # get sample paths
    s <- c(unlist(strsplit(list_sampleNames[j,2],split=',', fixed=TRUE)))
    samples <- lapply(s, function(x) if(which(s == x) %in% vector) return(x))
    samples <- plyr::compact(samples)
    samplePaths <- c(samplePaths,paste(data_directory, sample, samples, paste(samples,"_rawcounts.tab",sep=""), sep="/"))
    samplePaths_TPM <- c(samplePaths_TPM,paste(data_directory, sample, samples, paste(samples,"_tmp.tab",sep=""), sep="/"))
    
     # group conditions
    cond <- list_sampleNames[j,1]
    groupconditions <- c(groupconditions,rep(cond,length(samples)))
    
    # replicates
    s <- c(unlist(strsplit(list_sampleNames[j,4],split=',', fixed=TRUE)))
    repl <- lapply(s, function(x) if(which(s == x) %in% vector) return(as.numeric(x)))
    repl <- plyr::compact(repl)
    replicates <- c(replicates,unlist(repl))
  }
  
  ####################################### NOISeq object - raw FPKM + filtering ##############################################
  # save TPM
  mycounts <- data.frame()
  
  for(k in 1:length(samplePaths_TPM)){
    print(k)
    if (k==1){
      aux <- data.frame(read.csv2(file = samplePaths_TPM[[k]], sep="\t", header = T))
      mycounts <- data.frame(as.numeric(as.character(aux[,2])))
      colnames(mycounts) <- colnames(aux)[2]
      rownames(mycounts) <- aux[,1]
    }else{
      aux <- data.frame(read.csv2(file = samplePaths_TPM[[k]], sep="\t", header = T))
      m <- data.frame(as.numeric(as.character(aux[,2])))
      colnames(m) <- colnames(aux)[2]
      rownames(m) <- aux[,1]
      merging <- merge(mycounts,m,by.x = 0, by.y = 0,x.all=F,y.all=F)
      mycounts <- merging[,2:dim(merging)[2]]
      rownames(mycounts) <- merging$Row.names
    }
  } 
  filename <- paste(data_directory,sample,"tpm.tab",sep="/")
  write.table(mycounts,file=filename,row.names = T, col.names = T, quote = F,sep="\t")
  
  
  
  print("create noiseq obj ...")
  # my counts for NOISeq
  mycounts <- data.frame()
  
  for(k in 1:length(samplePaths)){
    print(k)
    if (k==1){
      aux <- data.frame(read.csv2(file = samplePaths[[k]], sep="\t", header = T))
      mycounts <- data.frame(as.numeric(as.character(aux[,2])))
      colnames(mycounts) <- colnames(aux)[2]
      rownames(mycounts) <- aux[,1]
    }else{
      aux <- data.frame(read.csv2(file = samplePaths[[k]], sep="\t", header = T))
      m <- data.frame(as.numeric(as.character(aux[,2])))
      colnames(m) <- colnames(aux)[2]
      rownames(m) <- aux[,1]
      merging <- merge(mycounts,m,by.x = 0, by.y = 0,x.all=F,y.all=F)
      mycounts <- merging[,2:dim(merging)[2]]
      rownames(mycounts) <- merging$Row.names
    }
  } 
  
  filename <- paste(data_directory,sample,"rawcounts.tab",sep="/")
  write.table(mycounts,file=filename,row.names = T, col.names = T, quote = F,sep="\t")
  
  # get design matrix
  myfactors <- data.frame(State=groupconditions, Replicates=replicates)
  row.names(myfactors) <- colnames(mycounts)
  
  ### format annotations
  mychroms <- data.frame(anno[,c(1,2,3)])
  row.names(mychroms) <- anno$V5
  colnames(mychroms) <- c("Chr","GeneStart","GeneEnd")
  
  ## biotypes
  mybiotypes <- substr(row.names(mychroms),1,3)
  names(mybiotypes) <- c(row.names(mychroms))
  
  

  ######### filtering ################3
  myfilt = filtered.data(mycounts, factor = myfactors$State, norm = FALSE, depth = NULL, method = 3,cpm = 5, p.adj = "fdr")
  mydata <- readData(data=myfilt,biotype=mybiotypes,chromosome=mychroms, factors=myfactors)
  
  myPCA = dat(mydata, type = "PCA")
  
  filename=paste(working_dir,"before_batch_removal.png",sep="/")
  png(file=filename,width=800,height=400)
  par(mfrow = c(1, 2))
  explo.plot(myPCA, factor = "State")
  explo.plot(myPCA, factor = "Replicates")
  dev.off()
  
  
  mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 3)
  explo.plot(mycd)
  
  ############# batch effect
  
  mydata2corr1 = ARSyNseq(mydata, factor = "Replicates", batch = TRUE, norm = "uqua", logtransf = FALSE)
  myPCA = dat(mydata2corr1, type = "PCA")
  
  filename=paste(working_dir,"after_batch_removal.png",sep="/")
  png(file=filename,width=800,height=400)
  par(mfrow = c(1, 2))
  explo.plot(myPCA, factor = "State")
  explo.plot(myPCA, factor = "Replicates")
  dev.off()
  
  mycd = dat(mydata2corr1, type = "cd", norm = FALSE, refColumn = 3)
  explo.plot(mycd)
  
  
  ############ remove noise
  par(mfrow = c(1, 2))
  mynoiseq1 = noiseqbio(mydata, k = 0.0001, norm = "uqua", plot=TRUE, factor="State",conditions=c("After","Before"),r=100)
  mynoiseq2 = noiseqbio(mydata2corr1, k = 0.0001, norm = "uqua", plot=TRUE, factor="State",conditions=c("After","Before"),r=100)
  
  results<- (mynoiseq2@results[[1]])
  mynoiseq2.deg = degenes(mynoiseq2, q = 0.95, M = NULL)
  mynoiseq2.deg1 = degenes(mynoiseq2, q = 0.95, M = "up")
  
  DE.plot(mynoiseq2, q = 0.95, graphic = "expr", log.scale = TRUE)
  DE.plot(mynoiseq2, q = 0.95, graphic = "MD")
  
  DE.plot(mynoiseq2, chromosomes = c(1, 2,16), log.scale = TRUE, join = FALSE, q = 0.85, graphic = "chrom")
  
  # myexplodata <- dat(mydata2corr1, type = "biodetection")
  
  # for(i in 1:subset_size){
  #     par(mfrow=c(1,2))
  #     explo.plot(myexplodata, samples=c(i,i+subset_size), plottype = "persample")
  # }
  # mynicedata <- dat2save(myexplodata)
  
  ####################################### NOISeq object - gene expression #######################################
  ##################################################################################################################
  
  name <- "geneexp_noiseq.tab"
  filename <- paste(working_dir,name,sep="/")
  
  write.table(mynoiseq2.deg,file=paste(working_dir,name,sep="/"),row.names = T, col.names = T, quote = F,sep="\t")
  
  
  ####################################### Ballgown object - transcript expression ##################################
  ##################################################################################################################
  
}
