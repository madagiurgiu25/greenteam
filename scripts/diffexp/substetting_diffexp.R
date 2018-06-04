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


