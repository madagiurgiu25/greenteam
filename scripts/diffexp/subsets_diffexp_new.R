setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown/Subsets")
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
set.seed(1234)

filterOne <- function(expr_data, na.rm = TRUE){
  size <- dim(expr_data)[2]
  results <- c(rep(TRUE,dim(expr_data)[1]))
  for (i in 1:dim(expr_data)[1]){
    for(j in 1:size)
      if(expr_data[i,j] < 1 ){
        results[i] <- FALSE
      }
  }
  return(results)
}

filterInfinite <- function(expr_data, na.rm = TRUE){
  size <- dim(expr_data)[2]
  results <- c(rep(TRUE,dim(expr_data)[1]))
  for (i in 1:dim(expr_data)[1]){
    for(j in 1:size)
      if(expr_data[i,j] > 10000 ){
        results[i] <- FALSE
      }
  }
  return(results)
}


filterHigh <- function(expr_data, na.rm = TRUE){
  size <- dim(expr_data)[2]
  results <- c(rep(TRUE,dim(expr_data)[1]))
  for (i in 1:dim(expr_data)[1]){
    for(j in 1:size)
      if(expr_data[i,j] > 500 ){
        results[i] <- FALSE
      }
  }
  return(results)
}


filterZero <- function(expr_data, na.rm = TRUE){
  size <- dim(expr_data)[2]
  results <- c(rep(TRUE,dim(expr_data)[1]))
  for (i in 1:dim(expr_data)[1]){
    for(j in 1:size)
      if(expr_data[i,j] ==0 ){
        results[i] <- FALSE
      }
  }
  return(results)
}

filterZero2 <- function(expr_data, na.rm = TRUE){
  size <- dim(expr_data)[2]
  results <- c(rep(TRUE,dim(expr_data)[1]))
  for (i in 1:dim(expr_data)[1]){
    count <- 0
    for(j in 1:(size/2))
      if(expr_data[i,j] !=0 ){
        count <- count +1
      }
    if (count< (size/2 -1)){
      results[i] <- FALSE
    }else{
      count <- 0
      for(j in (size/2+1):size)
        if(expr_data[i,j] !=0 ){
          count <- count +1
        }
        if (count< (size/2 -1)){
          results[i] <- FALSE
        }
    }
  }
  return(results)
}

filterFPKM <- function(expr_data, na.rm = TRUE){
  size <- dim(expr_data)[2]
  results <- c(rep(TRUE,dim(expr_data)[1]))
  for (i in 1:dim(expr_data)[1]){
    for(j in 1:(size/2))
      if(expr_data[i,j] < expr_data[i,size/2 - j + 2] ){
        results[i] <- FALSE
      }
  }
  return(results)
}

filterFPKM_mean <- function(expr_data, na.rm = TRUE){
  size <- dim(expr_data)[2]
  results <- c(rep(TRUE,dim(expr_data)[1]))
  for (i in 1:dim(expr_data)[1]){
    half = size/2
    after = sum(expr_data[i,1:half])/half
    
    before = sum(expr_data[i,(half+1):size])/half
    if (before > after){
      results[i] <- FALSE
    }
    
  }
  return(results)
}

# input example (df, c(M103A1, M103A2), c(M103B1,M103B2), 2, A12, B12)
plotFPKMs <- function(df,after,before,subset_size,name_after,name_before){
    
    require(cowplot)
    theme_set(theme_cowplot(font_size=12)) # reduce default font size
    
    df_2 <- gather(df[,1:(2*size)],"Sample","FPKM",1:(2*size))
    df_2$Type <- substring(df_2$Sample,1,10)
    df_2$Replicate <- substring(df_2$Sample,11,12)
    
    df3_quantile <- df_2[as.numeric(as.character(df_2$FPKM)) < quantile(as.numeric(as.character(df_2$FPKM)), 0.95), ]
    
    p1<- ggplot(df3_quantile,aes(as.numeric(as.character(FPKM)),fill=Type)) + 
        geom_density(alpha = 0.4) + 
         xlab("FPKM") + ggtitle(paste("FPKM distribution ",sample,sep=""))
    
    p2<-ggplot(df3_quantile,aes(x=Type, y=as.numeric(as.character(FPKM)),fill=Replicate)) + 
        geom_boxplot(alpha = 0.4) + scale_y_log10() + ylab("log10(FPKM)") + xlab("Condition") + 
        ggtitle(paste("FPKM distribution ",sample,sep=""))
    
    p3<-ggplot(df3_quantile, aes(x=as.numeric(as.character(FPKM)), fill=Replicate)) +
            geom_density(alpha=0.4) +
        facet_grid(~Type) + xlab("FPKM") +  ylab("density") + 
        ggtitle(paste("FPKM density ",sample,sep=""))
    
    png(file=paste("FPKM",name_after,paste(name_before,".png",sep=""),sep="_"),height=1200,width=800)
    grid.arrange(p1,p2,p3,ncol = 1)
    dev.off()
    
}


################################################### run subsets ###########################################################
############################# unmatched/matched ###########################################################################

data_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown')
data_directory_bam = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/STAR')

runBallgownSubsets <- function(sample,size, subset_size, unmatched){
  
  # change workingdir
  setwd(paste(data_directory,"Subsets",sample,sep="/"))
  
  sampleList<-combn(1:size, subset_size, simplify = FALSE)
  
  for (i in 1:length(sampleList)){
    
    vector <- unlist(sampleList[i], use.names = F)
    
    after <- paste(paste(sample, "A", sep=""),vector, sep="")
    after_samples<- paste("A",paste(vector, collapse = ''), sep="")
    
    before <- paste(paste(sample, "B", sep=""),vector, sep="")
    before_samples<- paste("B",paste(vector, collapse = ''), sep="")
    
    sample_IDs <- c(after,before)
    print(sample_IDs)
    
    # read samples into ballgown object
    sample_paths = paste(data_directory, sample_IDs, sep="/")
    print(sample_paths)
    bam_paths = paste(data_directory_bam, sample_IDs, "star_sorted.bam", sep="/")
    print(bam_paths)
    bg = ballgown(samples=sample_paths, bamfiles=bam_paths, meas='all')
    
    # matrix design
    pData(bg) = data.frame(id=sampleNames(bg), group=factor(rep(2:1, c(subset_size,subset_size))))
    print(pData(bg)) 
    
    # gene expression
    expression_genes = gexpr(bg)

    # all FPKM > 0
    fone=filterZero(expression_genes)
    df_genes_expression <- data.frame(expression_genes)
    df_genes_expression$filterZero=fone
    dim(df_genes_expression[df_genes_expression$filterZero== TRUE,])
    
    # plots
    plotFPKMs(subset(df_genes_expression[df_genes_expression$filterZero== TRUE,]), after, before, subset_size,after_samples,before_samples)
    
    fi=row.names(subset(df_genes_expression, df_genes_expression$filterZero == TRUE))
    bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)
    
    if (unmatched == TRUE){
      adjusted_results = stattest(bg_filtered, pData=pData(bg_filtered), feature='gene', meas='FPKM', covariate="group", getFC=TRUE)
    }else{
      matched_pair= factor(rep(1:subset_size,2))
      matched_pair
      mod = model.matrix(~matched_pair + pData(bg_filtered)$group)
      mod0 = model.matrix(~ pData(bg_filtered)$group)
      
      adjusted_results = stattest(bg_filtered, pData=pData(bg_filtered), feature='gene', meas='FPKM', mod0=mod0, mod=mod)
    }
    
    results_genes <- arrange(adjusted_results, qval)
    
    if (unmatched == TRUE){
      #filter only lnc significant
      re_lnc <- subset(results_genes, results_genes$qval < 0.05 & substring(results_genes$id,1,3) == 'LNC')
      expression_lnc <- subset(df_genes_expression, row.names(df_genes_expression) %in% re_lnc$id)
      write.table(expression_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_fpkm.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = T)
      write.table(re_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_logFC.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = F)
      
      #filter enriched in After IP
      expression_lnc <- subset(df_genes_expression, df_genes_expression$filterFPKM_mean == TRUE & row.names(df_genes_expression) %in% re_lnc$id)
      re_lnc <- subset(re_lnc,re_lnc$id %in% row.names(expression_lnc))
      write.table(expression_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_fpkm_enrichedAfterIP.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = T)
      write.table(re_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_logFC_enrichedAfterIP.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = F)
    }else{
      #filter only lnc significant
      re_lnc <- subset(results_genes, results_genes$qval < 0.05 & substring(results_genes$id,1,3) == 'LNC')
      expression_lnc <- subset(df_genes_expression, row.names(df_genes_expression) %in% re_lnc$id)
      write.table(expression_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_fpkm_matched.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = T)
      write.table(re_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_logFC_matched.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = F)
      
      #filter enriched in After IP
      expression_lnc <- subset(df_genes_expression, df_genes_expression$filterFPKM_mean == TRUE & row.names(df_genes_expression) %in% re_lnc$id)
      re_lnc <- subset(re_lnc,re_lnc$id %in% row.names(expression_lnc))
      write.table(expression_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_fpkm_enrichedAfterIP_matched.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = T)
      write.table(re_lnc,file=paste(sample,after_samples,before_samples,"expressed_genes_logFC_enrichedAfterIP_matched.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names = F)
      
    }
  }
  
}

############################## M103
# run M103 with subsamples of 4 + unmatched
runBallgownSubsets("M103",5,4,TRUE)

# run M103 with subsamples of 4 + matched
runBallgownSubsets("M103",5,4,FALSE)

# run M103 with subsamples of 3 + unmatched
runBallgownSubsets("M103",5,3,TRUE)

# run M103 with subsamples of 3 + matched
runBallgownSubsets("M103",5,3,FALSE)

# run M103 with subsamples of 5 + unmatched
runBallgownSubsets("M103",5,5,TRUE)

# run M103 with subsamples of 5 + matched
runBallgownSubsets("M103",5,5,FALSE)

############################# Mlet7

# run M103 with subsamples of 5 + unmatched
runBallgownSubsets("Mlet7",7,5,TRUE)

# run M103 with subsamples of 5 + unmatched
runBallgownSubsets("Mlet7",7,5,FALSE)


################################### Plot FPKMs
library(tidyr)
setwd("E:\\masterpraktikum\\mouse")

# df <- subset(df_genes_expression[df_genes_expression$filterZero== TRUE,])
df<-read.csv2(file="test.txt",sep = "\t",row.names=1,header=T)
size=4
df_2 <- gather(df[,1:8],"Sample","FPKM",1:(2*size))
df_2$Type <- substring(df_2$Sample,1,10)
df_2$Replicate <- substring(df_2$Sample,11,12)

ggplot(df_2,aes(as.numeric(as.character(FPKM)),fill=Type)) + 
    geom_density(alpha = 0.4) + scale_x_continuous( trans='log2') + 
    xlim(0,50) + scale_fill_manual(values=c("After","Before")) +
    xlab("FPKM") + ggtitle(paste("FPKM distribution ",sample,sep=""))

ggplot(df_2,aes(x=Type, y=as.numeric(as.character(FPKM)),fill=Replicate)) + 
    geom_boxplot(alpha = 0.4) + scale_y_log10() + ylab("log10(FPKM)") + xlab("Condition") + 
     ggtitle(paste("FPKM distribution ",sample,sep=""))

df<-read.csv2(file="overlaps_lnc.txt",sep="\t",header=F)
ggplot(df, aes(x=as.numeric(as.character(V1)))) + geom_step(aes(y=..y..),stat="ecdf",size=1) +
    xlab("Overlaps") + ylab("cummulative frequency of overlaps") + ggtitle("Cummulative plot LNC - overlaps") + scale_x_log10()


