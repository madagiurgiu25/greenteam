setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown")
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

########################### new filtering + clustering
sample="M103"
size=5

# read samples into R
data_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown')
sample_IDs = list.dirs(path = data_directory, full.names = FALSE, recursive = FALSE)
sample_IDs <- sample_IDs[ grepl(paste("^",sample,sep=""), sample_IDs) ]
print(sample_IDs)
sample_paths = paste(data_directory, sample_IDs, sep="/")


# bam files
data_directory_bam = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/STAR')
bam_paths = paste(data_directory_bam, sample_IDs, "star_sorted.bam", sep="/")


# create ballgown object
bg = ballgown(samples=sample_paths, bamfiles=bam_paths, meas='all')

# get annotation
#anno_path = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/mapping/mm10_primary_assembly_and_lncRNA.gtf')
#annot = gffReadGR(anno_path, splitByTranscript=TRUE)
#info = annotate_assembly(assembled=structure(bg)$trans, annotated=annot)

# design matrix
pData(bg) = data.frame(id=sampleNames(bg), group=factor(rep(2:1, c(size,size))))
pData(bg) 
# mark the matched paired
matched_pair= factor(rep(1:size,2))
matched_pair
mod = model.matrix(~matched_pair + pData(bg)$group)
mod0 = model.matrix(~ pData(bg)$group)

# get GENE EXPRESSION level
expression_genes = gexpr(bg)
l<-cbind(t_id=texpr(bg, 'all')$t_name,gene_id=texpr(bg, 'all')$gene_id)
write.table(l,file="count_transcripts_per_gene.txt",sep="\t",quote=F,row.names=F)

#filter my gene expression
df_genes_expression = data.frame(expression_genes)

# all FPKM > 1
fone=filterOne(expression_genes)
df_genes_expression$filterOne=fone
dim(df_genes_expression[df_genes_expression$filterOne == TRUE,])

# all FPKM > 0
fone=filterZero(expression_genes)
df_genes_expression$filterZero=fone
dim(df_genes_expression[df_genes_expression$filterZero== TRUE,])

map_transcript2gene<-data.frame(t_id=texpr(bg, 'all')$t_name,gene_id=texpr(bg, 'all')$gene_id)
count_t_per_gene<-(ddply(map_transcript2gene,.(gene_id),nrow))
df_transcript2_per_gene <- data.frame(frequency=as.numeric(as.character(count_t_per_gene[,2])), gene=count_t_per_gene[,1], type=(substring(count_t_per_gene[,1],1,3)))


tr <- c(100,10,50,8)
k <- c(20,3,10,3)
i=1
for (i in 1:5){
    
    print("read data")
    # cluster expression
    tr_aux <- tr[i]
    k_aux <- k[i]
    df <- read.table(file=paste("Clustering/",sample,"_cluster_tr_",tr_aux,"_k_", k_aux, "_.txt",sep=""),header=T, row.names=1,sep="\t")
    
    print("append genes which are not clustered")
    # the rest of the genes
    df_genes=subset(df_transcript2_per_gene, df_transcript2_per_gene$frequency <= tr_aux )
    list_genes=as.character(df_genes$gene)
    bg_filter = subset(bg,"ballgown::geneIDs(bg) %in% list_genes",genomesubset=TRUE)
    expression_matrix = gexpr(bg_filter)
    # row.names(expression_matrix) <- ballgown::geneIDs(bg_filter)
    df <- rbind(df, expression_matrix)  
    
    print("stattest FPKM >1 ")
    # statistical test FPKM > 1
    df_filtered = subset(df, filterOne(df) == TRUE & filterInfinite(df) == TRUE)
    
    adjusted_results = stattest(gowntable = df_filtered, pData=pData(bg), feature='gene', meas='FPKM', mod0=mod0, mod=mod)
    results_transcript <- arrange(adjusted_results, qval)
    print(head(results_transcript,50))
    
    print("stattest FPKM >0 ")
    # statistical test FPKM > 0
    df_filtered = subset(df, filterZero(df) == TRUE & filterInfinite(df) == TRUE)
    
    adjusted_results = stattest(gowntable = df_filtered, pData=pData(bg), feature='gene', meas='FPKM', mod0=mod0, mod=mod)
    results_transcript <- arrange(adjusted_results, qval)
    print(head(results_transcript,50))
    }

#############################################################################
###################################### group only some replicate - not paired
data_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown')
sample_IDs <- c("M103A1",
                "M103A2",
                "M103A3",
                "M103A5",
                "M103B1",
                "M103B2",
                "M103B5",
                "M103B4"
)
print(sample_IDs)
sample_paths = paste(data_directory, sample_IDs, sep="/")


# bam files
data_directory_bam = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/STAR')
bam_paths = paste(data_directory_bam, sample_IDs, "star_sorted.bam", sep="/")


# create ballgown object
bg = ballgown(samples=sample_paths, bamfiles=bam_paths, meas='all')

# mark the matched paired

expression_genes = gexpr(bg)
# filter my gene expression
df_genes_expression = data.frame(expression_genes)

# all FPKM > 1
fone=filterOne(expression_genes)
df_genes_expression$filterOne=fone
dim(df_genes_expression[df_genes_expression$filterOne == TRUE,])


# all FPKM > 0
fzero=filterZero(expression_genes)
df_genes_expression$filterZero=fzero
dim(df_genes_expression[df_genes_expression$filterZero == TRUE,])

# FPKM_a > FPKM_b
fzero=filterFPKM(expression_genes)
df_genes_expression$filterFPKM=fzero
dim(df_genes_expression[df_genes_expression$filterFPKM == TRUE,])
head(df_genes_expression[df_genes_expression$filterFPKM == TRUE,])

# mean FPKM_a > FPKM_b
fzero=filterFPKM_mean(expression_genes)
df_genes_expression$filterFPKM_mean=fzero
dim(df_genes_expression[df_genes_expression$filterFPKM_mean == TRUE,])

# mean FPKM <= 30000
fzero=filterInfinite(expression_genes)
df_genes_expression$filterInfinite =fzero
dim(df_genes_expression[df_genes_expression$filterInfinite  == TRUE,])



fi=row.names(subset(df_genes_expression, df_genes_expression$filterZero == TRUE))

pData(bg) = data.frame(id=sampleNames(bg), group=factor(rep(2:1, c(4,4))))
pData(bg) 
bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

adjusted_results = stattest(bg_filtered, pData=pData(bg_filtered), feature='gene', meas='FPKM', covariate="group", getFC=TRUE)
results_genes <- arrange(adjusted_results, qval)

#filter only lnc significant
re_lnc <- subset(results_genes, results_genes$qval < 0.05 & substring(results_genes$id,1,3) == 'LNC')
#get expression and filter only those lnc mean(FPKM_B) < mean(FPKM_A)
expression_lnc <- subset(df_genes_expression, df_genes_expression$filterFPKM_mean == TRUE & row.names(df_genes_expression) %in% re_lnc$id)
write.table(expression_lnc,file="M103_A1235_B1245_expressed_genes_significant.txt",sep="\t",col.names=T,quote=F,row.names = F)
write.table(results_genes,file="M103_A1235_B1245_expressed_genes.txt",sep="\t",col.names=T,quote=F,row.names = F)


##########################################################################################################################
################################################## Plot replicates #######################################################
library(reshape)

log2FC <- function(col_before, col_after){
  size <- length(col_before)
  results <- c(rep(0,size))
  for (i in 1:size){
    if (col_after[i] == col_before[i] && col_after[i] == 0){
      results[i] <- 0
    }else if (col_before[i] == 0){
      results[i] <- log2((col_after[i]/(col_before[i]+1)))
    }else{
      results[i] <- log2(col_after[i]/(col_before[i]))
    }
  }
  return(results)
}

checkFC <- function(df, treshold){
  size <- dim(df)[2]
  length <- dim(df)[1]
  results <- rep(0,length)
  for (i in 1:length){
    count = 0
    for (j in 1:size){
      if (abs(df[i,j]) > treshold){
        count <- count + 1
      }
    results[i] <- count
    }
  }
  return(results)
}


sample="M103"
size=5
data_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown')
sample_IDs <- c("M103A1",
                "M103A2",
                "M103A3",
                "M103A4",
                "M103A5",
                "M103B1",
                "M103B2",
                "M103B3",
                "M103B4",
                "M103B5"
)
sample_paths = paste(data_directory, sample_IDs, sep="/")
# bam files
data_directory_bam = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/STAR')
bam_paths = paste(data_directory_bam, sample_IDs, "star_sorted.bam", sep="/")

# create ballgown object
bg = ballgown(samples=sample_paths, bamfiles=bam_paths, meas='all')

# mark the matched paired

expression_genes = gexpr(bg)
# filter my gene expression
df_genes_expression = data.frame(expression_genes)
fzero=filterZero(expression_genes)
df_genes_expression$filterZero=fzero

## genes with FPKM > 0
fi=row.names(subset(df_genes_expression, df_genes_expression$filterZero == TRUE))

pData(bg) = data.frame(id=sampleNames(bg), group=factor(rep(2:1, c(size,size))))
pData(bg) 

### filtered diff expr
filtered_df <- subset(df_genes_expression, df_genes_expression$filterZero == TRUE)
filtered_df$id <- row.names(filtered_df)

after <- colnames(filtered_df)[1:size]
before <- colnames(filtered_df)[size+1:(2*size)]

# after <- c("FPKM.M103A1" ,"FPKM.M103A2" ,"FPKM.M103A4" ,"FPKM.M103A3", "FPKM.M103A5")
# before <- c("FPKM.M103B1" ,"FPKM.M103B2" ,"FPKM.M103B3" ,"FPKM.M103B4", "FPKM.M103B5")

for (i in 1:size){
  name <- paste(after[i],before[i],"log2FC",sep="_")
  filtered_df[[name]] <- log2FC(filtered_df[[before[i]]],filtered_df[[after[i]]])
}

filtered_df$checklogFC <- checkFC(filtered_df[,13:17],2)

bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)
adjusted_results = stattest(bg_filtered, feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
results_genes = arrange(adjusted_results, qval)


########### plots
library(reshape)
df_logfc <- filtered_df[,13:17]
mdata <- melt(df_logfc) 

ggplot(mdata, aes(x=as.numeric(as.character(value)),color=variable)) + geom_step(aes(y=..y..),stat="ecdf") +
  xlab("log2FC") + ylab("cummulative frequency of log2FC") + ggtitle("Cummulative plot M103 log2FC (-5,5)") + 
  xlim(-5,5)

ggplot(mdata, aes(x=as.numeric(as.character(value)),color=variable)) + geom_step(aes(y=..y..),stat="ecdf") +
  xlab("log2FC") + ylab("cummulative frequency of log2FC") + ggtitle("Cummulative plot M103 log2FC") 


df_plot <- data.frame(filtered_df, type=substring(row.names(filtered_df),1,3))
ggplot(data=df_plot, aes(x=as.character(checklogFC),fill=type)) +
  geom_bar(position=position_dodge()) +
  geom_text(stat='count', aes(label=..count..), position = position_dodge(width = 1),
            vjust = -0.5) + 
  xlab("#replicates with |log2FC| > 2") + ggtitle("M103 expression consistency over genes")


########################################################## add logFC from ballgown

filtered_df_merged <- merge(x=filtered_df,y=results_genes, by="id")
filtered_df_merged_sort<- filtered_df_merged[order(filtered_df_merged$qval),]

setwd('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown/Subsets/M103')
row.names(filtered_df_merged_sort) <- filtered_df_merged_sort$id
write.table(filtered_df_merged_sort,file="logFC_FPKM_qval_M103.txt",sep="\t",col.names = T,quote=F,row.names = F )

lnc_group_345<-subset(filtered_df_merged_sort, filtered_df_merged_sort$checklogFC >=3 & substring(row.names(filtered_df_merged_sort),1,3) == 'LNC')
write.table(lnc_group_345,file="logFC_FPKM_qval_M103_onlyLNC_class345.txt",sep="\t",col.names = T,quote=F,row.names = F )

lnc_group_45<-subset(filtered_df_merged_sort, filtered_df_merged_sort$checklogFC >3 & substring(row.names(filtered_df_merged_sort),1,3) == 'LNC')
write.table(lnc_group_45,file="logFC_FPKM_qval_M103_onlyLNC_class45.txt",sep="\t",col.names = T,quote=F,row.names = F )

lnc_group_45<-subset(filtered_df_merged_sort, filtered_df_merged_sort$checklogFC >3)
write.table(lnc_group_45,file="logFC_FPKM_qval_M103_class45.txt",sep="\t",col.names = T,quote=F,row.names = F )

lnc_group_45<-subset(filtered_df_merged_sort, filtered_df_merged_sort$checklogFC >=3)
write.table(lnc_group_45,file="logFC_FPKM_qval_M103_class345.txt",sep="\t",col.names = T,quote=F,row.names = F )



###### boxplots
df_plot <- data.frame(filtered_df_merged_sort, type=substring(row.names(filtered_df_merged_sort),1,3))
ggplot(df_plot, aes(x = as.factor(checklogFC), y = -log10(as.numeric(as.character(qval))),fill=type)) +
  geom_boxplot() + ylab("-log10(qval)") + xlab("#count replicates with |log2FC| > 2") + ggtitle("qval distribution for |logFC|>2")

########## ballgown fc + replicates pairs logFC
df_logfc <- filtered_df_merged_sort[,c(13,14,15,16,17,20)]
df_logfc$fc <- log2(df_logfc$fc)
mdata <- melt(df_logfc) 
ggplot(mdata, aes(x=as.numeric(as.character(value)),color=variable)) + geom_step(aes(y=..y..),stat="ecdf",size=1) +
  # scale_color_manual(values=c("green", "magenta", "#56B4E9","red","blue","black"),name="replicates") + 
  scale_color_manual(labels = c("log2FC(A1-B1)", "log2FC(A2-B2)","log2FC(A3-B3)","log2FC(A4-B4)","log2FC(A5-B5)","log2FC_ballgown"),values=c("green", "magenta", "#56B4E9","red","blue","black"),name="replicates") + 
  xlab("log2FC") + ylab("cummulative frequency of log2FC") + ggtitle("Cummulative plot M103 log2FC") + xlim(-10,10)


######################## ONLY LNC
df_logfc <- filtered_df_merged_sort[,c(13,14,15,16,17,20)]
df_logfc$fc <- log2(df_logfc$fc)
df_lnc <- subset(df_logfc, substring(row.names(df_logfc),1,3) == 'LNC')
mdata <- melt(df_lnc) 
ggplot(mdata, aes(x=as.numeric(as.character(value)),color=variable)) + geom_step(aes(y=..y..),stat="ecdf",size=1) +
  scale_color_manual(labels = c("log2FC(A1-B1)", "log2FC(A2-B2)","log2FC(A3-B3)","log2FC(A4-B4)","log2FC(A5-B5)","log2FC_ballgown"),values=c("green", "magenta", "#56B4E9","red","blue","black"),name="replicates") + 
  xlab("log2FC") + ylab("cummulative frequency of log2FC") + ggtitle("Cummulative plot M103 log2FC (lncRNA only)") + xlim(-10,10)


################### FPKM cummulative
df_logfc <- filtered_df_merged_sort[,c(2:6)]
mdata <- melt(df_logfc) 
ggplot(mdata, aes(x=as.numeric(as.character(value)),color=variable)) + geom_step(aes(y=..y..),stat="ecdf",size=1) +
  scale_color_manual(labels = c("A1", "A2","A3","A4","A5"),values=c("green", "magenta", "#56B4E9","red","blue"),name="samples After IP") + 
  xlab("log10(FPKM)") + ylab("cummulative frequency of FPKM") + ggtitle("Cummulative plot M103 FPKM") + scale_x_log10()

# tab<-data.frame(table(mdata$value))
# tab = transform(tab,cumsum=cumsum(Freq))
# tag_merge<-merge(tab,mdata, by.x="Var1",by.y="value")
# ggplot(tag_merge,aes(y=as.numeric(as.character(cumsum)),x=as.numeric(as.character(Var1)),color=variable)) + geom_line() + scale_x_log10()


df_logfc <- filtered_df_merged_sort[,c(7:11)]
mdata <- melt(df_logfc) 
ggplot(mdata, aes(x=as.numeric(as.character(value)),color=variable)) + geom_step(aes(y=..y..),stat="ecdf",size=1) +
  scale_color_manual(labels = c("B1", "B2","B3","B4","B5"),values=c("green", "magenta", "#56B4E9","red","blue"),name="samples Before IP") + 
  xlab("log10(FPKM)") + ylab("cummulative frequency of FPKM") + ggtitle("Cummulative plot M103 FPKM") + scale_x_log10()


################################################# Plot Gene and lncRNA overlaps FC expression
setwd('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown/Subsets/M103/consistent_lncRNA_class45')
data <- read.csv2(file="lnc_ens_overlaps_class45_enrichedAfterIP.txt",sep="\t", header=F)
colnames(data)<-c("lnc","strand","gene","strand","overlap_length","lnc_fraction","gene_fraction","type")

data_short<-data.frame(lnc=data$lnc, gene=data$gene,lnc_fraction=data$lnc_fraction)

# get all ballgown logFC for lncRNA
aux <- subset(filtered_df_merged_sort, filtered_df_merged_sort$id %in% data$lnc)
df_lnc <- data.frame(lnc=aux$id,fc_lnc=log2(aux$fc)) 
merge1 <- merge(data_short,df_lnc,by='lnc',x.all=TRUE)

# get all ballgown logFC for gene
aux <- subset(filtered_df_merged_sort, filtered_df_merged_sort$id %in% data$gene)
df_gene <- data.frame(gene=aux$id,fc_gene=log2(aux$fc),found=1) 

# genes which are not expressed or low expressed
genes_not_expressed <-  setdiff(data$gene,aux$id)
expr<- data.frame(subset(expression_genes, row.names(expression_genes) %in% genes_not_expressed))
df_low_genes<- data.frame(gene=row.names(expr), fc_gene = log2FC(rowMeans(expr[,1:5]),rowMeans(expr[,6:10])),found=0)

df_gene<- rbind(df_gene,df_low_genes)

merge2 <- merge(merge1,df_gene,by='gene',x.all=TRUE)

ggplot(merge2, aes(x=fc_lnc,y=fc_gene, shape=as.factor(found), color=as.numeric(as.character(lnc_fraction)))) + geom_point() + 
  xlab("lnc log2FC") + ylab("gene log2FC") +
  scale_colour_gradient(name="lnc overlapping fraction",low = "red", high = "gray") +
  scale_shape_discrete(labels=c("low-expressed","diff-exp"),name="gene types") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + ggtitle("logFC correlation between lncRNA and its overlapping gene \n (lnc = consistent diff exp over at least 4 replicates and enriched in After IP)")


df_subset<-subset(merge2,merge2$found == 1 & merge2$fc_lnc >= 2 & abs(merge2$fc_gene)<1)
write.table(df_subset,file="subset_lnc_gene_overlap_fclnc>2_fcgene<1_both_expressed.txt",sep="\t",col.names = T,quote=F,row.names = F)


