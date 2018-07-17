#!/usr/bin/env Rscript
library(ballgown)
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
# library(here)
library(dplyr)


setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/scripts/diff_exp/")
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
data_bam_directory = (args[4])
working_dir = (args[5])
printGenes = as.logical(args[6])
plots = as.logical(args[7])
transcriptexpression = as.logical(args[8])
# #samples condition
sampleNames_conditions = data.frame(read.csv2(file = args[9],sep = "\t",header = F))

# #genes with interassant meaning (2xcolumns {Gene_id   Type})
special_genes = data.frame(read.csv2(file = args[10],sep = "\t",header = F))


print(sample)
print(subset_size)
print(data_directory)
print(data_bam_directory)
print(working_dir)
print(printGenes)
print(plots)
print(transcriptexpression)
print(args[9])
print(args[10])

# 
# sample = "H103"
# subset_size = 6
# data_directory = (
# '/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown/H103'
# )
# data_bam_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/STAR/H103')
# printGenes = TRUE
# plots = TRUE
# transcriptexpression = TRUE
# working_dir = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/DiffExp/H103/')
# special_genes = data.frame(read.csv2(file = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/TextMining/human_mir103_formatted.txt_mapped" ,sep = "\t",header = F))
# # special_genes = data.frame(read.csv2(file = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/miranda/mm10/tmp.txt",sep = "\t",header = F))
# # 
# sampleNames_conditions = data.frame(
# read.csv2(file = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/samples/h103.txt", sep = "\t", header = F)
# )
# size = 6

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
  bamPaths <- vector()
  groupconditions <- vector()
  replicates <- vector()
  
  bamPathsString <- vector()
  
  for (j in 1:dim(list_sampleNames)[1]){
    
    ########################################## Format input #####################################################
    
    # get sample paths
    s <- c(unlist(strsplit(list_sampleNames[j,2],split=',', fixed=TRUE)))
    samples <- lapply(s, function(x) if(which(s == x) %in% vector) return(x))
    samples <- plyr::compact(samples)
    samplePaths <- c(samplePaths,paste(data_directory, samples, sep="/"))
    
    # get bam paths
    newbams <- paste(data_bam_directory, samples, "star_sorted.bam", sep="/")
    bamPaths <- c(bamPaths,newbams)
    bamPathsString <- c(bamPathsString,paste(newbams,collapse = ','))
    
    # group conditions
    cond <- list_sampleNames[j,1]
    groupconditions <- c(groupconditions,rep(cond,length(samples)))
    
    # replicates
    s <- c(unlist(strsplit(list_sampleNames[j,4],split=',', fixed=TRUE)))
    repl <- lapply(s, function(x) if(which(s == x) %in% vector) return(as.numeric(x)))
    repl <- plyr::compact(repl)
    replicates <- c(replicates,unlist(repl))
  }
  
  ####################################### Ballgown object - raw FPKM ##############################################
  print("create ballgown obj ...")
  # create ballgown object
  # bg_obj<-createBallgownObj(samplePaths,bamPaths)
  
  bg_obj = ballgown(samples = samplePaths, bamfiles=bamPaths, meas = 'all')
  print(head(ballgown::geneIDs(bg_obj)))
  
  ####################################### Ballgown object - gene expression #######################################
  ##################################################################################################################
  print("gene expression obj ...")
  # gene expression
  expression_genes = gexpr(bg_obj)
  
  #######################################        Filtering           ##############################################
  # all FPKM > 0
  print("filter FPKM > 0....")
  fone = filterZero(expression_genes)
  df_genes_expression <- data.frame(expression_genes)
  df_genes_expression$filterZero = fone
  fmean = filterFPKM_mean(expression_genes)
  df_genes_expression$filterFPKM_mean = fmean
  
  df_genes_expression <- subset(df_genes_expression, df_genes_expression$filterZero == TRUE)
  print(head(df_genes_expression))
  
  # store fpkm for further use
  if(printGenes){
    
    # write in file the FPKMs
    
    # path <- paste(working_dir,sample,sep="/")
    # name <- paste("fpkm",paste(replicates,collapse=''),sep="_")
    name <- "fpkm.tab"
    filename <- paste(working_dir,name,sep="/")
    print(paste("write gene expression to file: ",filename,sep=""))
    
    write.table(df_genes_expression[,1:length(replicates)],file=filename,row.names = T, col.names = T, quote = F,sep="\t")
  
    name <- "bams.txt"
    filename <- paste(working_dir,name,sep="/")
    write.table(paste(bamPathsString,collapse = ' '),file=filename,row.names = F,col.names = F,quote = F)
  }   
  
  if(plots){
    
    png(filename = paste(working_dir,"before_norm_hist.png",sep="/"),width = 800, height = 800)
        
    # plot distribution FPKM
    plotFPKM <- plotHistogram(df_genes_expression,subset_size,sample,groupconditions, replicates, "filtered")
    grid.arrange(plotFPKM)
    dev.off()
    
  }
  
  #######################################        Normalization           ##########################################
  
  ###Translate + Transform log + Quantile norm
  df_trans <- translate_transformlog2(df_genes_expression[,1:length(replicates)],subset_size,conditions)
  df_norm_1 <- quantile_normalization_perCondition_new(df_trans[,1:subset_size],subset_size)
  for (j in 2:conditions){
    df_norm_1 <- data.frame(df_norm_1, quantile_normalization_perCondition_new(df_trans[,((j-1)*(subset_size)+1):(j*subset_size)],subset_size))
  }

    ########### SQRT + Quantile Normalization
  df_trans <- transformSQRT(df_genes_expression[,1:length(replicates)],subset_size,conditions)
  df_norm_2 <- quantile_normalization_perCondition_new(df_trans[,1:subset_size],subset_size)
  for (j in 2:conditions){
    df_norm_2 <- data.frame(df_norm_2, quantile_normalization_perCondition_new(df_trans[,((j-1)*(subset_size)+1):(j*subset_size)],subset_size))
  }

  ########### Variance stablization + Quantile Normalization
  df_trans <- variance_stabilization(df_genes_expression[,1:length(replicates)],subset_size,conditions)
  df_norm_3 <- quantile_normalization_perCondition_new(df_trans[,1:subset_size],subset_size)
  for (j in 2:conditions){
    df_norm_3 <- data.frame(df_norm_3, quantile_normalization_perCondition_new(df_trans[,((j-1)*(subset_size)+1):(j*subset_size)],subset_size))
  }
  print("Normalization.....")
  print(head(df_norm_3))

  if(plots){
    
    png(filename = paste(working_dir,"after_norm_hist.png",sep="/"),width = 1600, height = 500)
    # plot distribution FPKM
    p1 <- plotHistogram(df_norm_1,subset_size,sample,groupconditions, replicates, "normalized 1 (log + quantile norm)")
    p2 <- plotHistogram(df_norm_2,subset_size,sample,groupconditions, replicates, "normalized 2 (sqrt + quantile norm)")
    p3 <- plotHistogram(df_norm_3,subset_size,sample,groupconditions, replicates, "normalized 3 (variance stabilization + quantile norm)")
    grid.arrange(p1,p2,p3, ncol=3)
    dev.off()
    
  }
  
  #######################################        Call ballgown            ##########################################
  df_norm_expr <- transformExp2(df_norm_3)
  # df_norm_expr <- df_norm_3
  print("Exp 2 .....")
  #all<-callBallgown_geneexpression(bg_obj, df_norm_3, subset_size, groupconditions,replicates)
  
  pData(bg_obj) = data.frame(id = sampleNames(bg_obj), group = groupconditions)
  print(pData(bg_obj))
  
  bg_filtered = subset(bg_obj, "ballgown::geneIDs(bg_obj) %in% row.names(df_norm_3)", genomesubset = TRUE)
  
  matrix <- as.matrix(mutate_all(df_norm_expr, function(x) as.numeric(as.character(x))))
  print(head(matrix))
  colnames(matrix) <- colnames(df_norm_expr)
  rownames(matrix) <- row.names(df_norm_expr)
  
  adjusted_results = stattest(
    gowntable = matrix,
    pData=pData(bg_filtered),
    feature = 'gene',
    meas = 'FPKM',
    covariate="group",
    getFC = TRUE,
    libadjust = FALSE
  )
  all <- arrange(adjusted_results, qval)
  
  print(head(all))
  subset_textmining<-data.frame(Names=special_genes$V1)
  subset_lnc_textmining<-subset(subset_textmining,substr(subset_textmining$Names,1,3) == 'LNC')
  
  # print vulcano
  if(plots){
    
    png(filename = paste(working_dir,"volcano.png",sep="/"),width = 800, height = 700)
    
    all$log2FoldChange<-log2(as.numeric(as.character(all$fc)))
    all$log10qval <- -log10(as.numeric(as.character(all$qval)))
    title<-paste(sample, "Diff.exp.genes", sep = " ")
    with(all, plot(log2FoldChange, log10qval, pch=20, main=title,ylab = "-log10(qval)",col="grey"))
    # with(subset(all, substring(all$id,1,3) == "LNC" ), points(log2FoldChange, log10qval, pch=20, col="black"))
    with(subset(all, as.numeric(as.character(all$qval)) < 0.05 & substring(all$id,1,3) == "LNC"), points(log2FoldChange, log10qval, pch=20, col="blue"))
    with(subset(all, substr(all$id,1,18) %in% subset_textmining$Names | all$id %in% subset_textmining$Name), points(log2FoldChange, log10qval, pch=20, col="red"))
    with(subset(all, all$id %in% subset_lnc_textmining$Names), points(log2FoldChange, log10qval, pch=20, col="green"))
    with(all,legend("topright",c("lncRNA, qval < 0.05", as.character(special_genes$V2[1]),paste("lncRNA",as.character(special_genes$V2[1]),sep=" ")), fill = c("blue","red","green")))
    dev.off()
    
  }
  genexp_name <- "genexp"
  if (size != subset_size){
    genexp_name <- paste(genexp_name, paste(replicates,collapse = ""), ".tab",sep = "")
    genexp_name_sign <- paste(genexp_name, paste(replicates,collapse = ""), "_signif.tab",sep = "")
  }else{
    genexp_name <- paste(genexp_name, ".tab",sep = "")
    genexp_name_sign <- paste(genexp_name, "_signif.tab",sep = "")
  }
  
  write.table(all,file=paste(working_dir,genexp_name,sep="/"),row.names = T, col.names = T, quote = F,sep="\t")
  write.table(subset(all, as.numeric(as.character(all$qval)) < 0.05),file=paste(working_dir,genexp_name_sign,sep="/"),row.names = T, col.names = T, quote = F,sep="\t")
  
  
  ####################################### Ballgown object - transcript expression ##################################
  ##################################################################################################################

}
