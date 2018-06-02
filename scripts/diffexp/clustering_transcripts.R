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
            if(expr_data[i,j] > 30000 ){
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
            if(expr_data[i,j] == 0 ){
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


sample=args[1]
size=as.integer(args[2])
transcript_per_gene_treshhold=as.integer(args[3])
k_clusters=as.integer(args[4])


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


#transcript and gene mapping
map_transcript2gene<-data.frame(t_id=texpr(bg, 'all')$t_name,gene_id=texpr(bg, 'all')$gene_id)
count_t_per_gene<-(ddply(map_transcript2gene,.(gene_id),nrow))
df_transcript2_per_gene <- data.frame(frequency=as.numeric(as.character(count_t_per_gene[,2])), gene=count_t_per_gene[,1], type=(substring(count_t_per_gene[,1],1,3)))



##############################################################################################################
##############################################################################################################
################################################### CLUSTERING ###############################################
##############################################################################################################
##############################################################################################################


# find genes that need to cluster
df_genes2cluster=subset(df_transcript2_per_gene, df_transcript2_per_gene$frequency > transcript_per_gene_treshhold )
list_genes2cluster=as.character(df_genes2cluster$gene)

# get empty format
df <- texpr(bg)
df <- df[df[,1] > 100000,]
mapping<-data.frame(cluster=integer(),
                    t_id=integer(), 
                    t_name=character(),
                    geneID=character(), 
                    stringsAsFactors=FALSE)

count = 0
# cluster genes
#for (i in 1:length(list_genes2cluster)){
for (i in 1:5){
    print(count)
    count = count + 1
    
    gene<-list_genes2cluster[i]
    print(gene)
    
    obj <- collapseTranscripts(gene, bg, meas = "FPKM", method = c("hclust"), k = k_clusters)
    
    #####3 add fpkm to matrix
    collapsed_transcripts <- obj$tab
    df_feature<- data.frame(cluster=c(paste("CLUST", gene, row.names(collapsed_transcripts), sep="_")))
    row.names(collapsed_transcripts) <- df_feature$cluster
    df<-rbind(df, collapsed_transcripts)
    
    
    ###### add cluster mapping
    transcripts<-data.frame(t_id=texpr(bg, 'all')$t_id, t_name=texpr(bg, 'all')$t_name)
    clust_assign <- df_feature[match(obj$cl$cluster$cluster,row.names(df_feature)),]
    clusters<- data.frame(t_id=obj$cl$clusters$t_id, geneID=rep(gene,dim(obj$cl$clusters)[1]),clusterName=clust_assign)
    # inner join
    merging<-merge(x=transcripts,y=clusters)
    mapping <- rbind(mapping,merging)
    
    
}
setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown/Clustering")
write.table(df,paste(sample,"cluster_tr", transcript_per_gene_treshhold ,"k", k_clusters, ".txt",sep="_"),sep="\t",row.names=T,quote=F)
write.table(mapping,paste(sample,"cluster_tr", transcript_per_gene_treshhold ,"k", k_clusters, "_mapping.txt",sep="_"),sep="\t",row.names=F,quote=F)
