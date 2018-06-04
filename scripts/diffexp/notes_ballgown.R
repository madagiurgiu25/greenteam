# read samples into R
data_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown')
sample_IDs = list.dirs(path = data_directory, full.names = FALSE, recursive = FALSE)
sample_IDs <- sample_IDs[ grepl("^M103", sample_IDs) ]
sample_paths = paste(data_directory, sample_IDs, sep="/")
bg = ballgown(samples=sample_paths, meas='all')

# access exon structure e_data
structure(bg)$exon
# access intron structure i_data
structure(bg)$intron
# access transcript structure t_data
structure(bg)$trans

plotTranscripts('NONMMUG034479.2', bg, samples=sample_IDs, meas='FPKM', colorby='transcript')

pData(bg) = data.frame(id=sampleNames(bg), group=rep(1:0, c(4,4))) # first 4 samples = 1, last 4 samples = 0
phenotype_table = pData(bg)
stat_results = stattest(bg, feature='gene', meas='FPKM', covariate='group', getFC=TRUE)
head(stat_results)

# Filter to remove low-abundance genes
bg_filtered = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)

# Identify transcripts that show statistically significant differences between groups.
results_transcripts = stattest(bg_filtered,feature="transcript", covariate="group", getFC=TRUE, meas="FPKM")

# Identify genes that show statistically significant differences between groups
results_genes = stattest(bg_filtered, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")

# Add gene names to results_transcripts
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filtered),
                                 geneIDs=ballgown::geneIDs(bg_filtered),
				 transcriptIDs=ballgown::transcriptNames(bg_filtered),
                                 results_transcripts)

# Sort the results from the smallest P value to the largest
results_transcripts = arrange(results_transcripts, pval)
results_genes = arrange(results_genes, pval)

head(results_genes, 20)
head(results_transcripts, 20)

tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)

# Show the distribution of gene abundances (measured as FPKM values) across samples, colored by group
fpkm = texpr(bg_filtered, meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm, col=as.numeric(phenotype_table$group), las=2, ylab='log2(FPKM+1)', main="Distribution of gene abundances across samples")

# Make plots of individual transcripts across samples
myTranscript <- ballgown::transcriptNames(bg_filtered)[2]
myGene <- ballgown::geneNames(bg_filtered)[2]
plot(fpkm[2,] ~ phenotype_table$group, border=c(1,2), main=paste(myGene,' : ',myTranscript),pch=19, xlab="Group", ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)),
        col=as.numeric(pheno_data$sex))

# Identify transcripts and genes with a q value <0.05
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

# View the distribution of differential expression values as a histogram. Display only those that are significant according to Ballgown
sig = which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
hist(results_genes[sig,"de"], breaks=50,  main="Distribution of differential expression values", col=as.numeric(phenotype_table$group), las=2)

# Plot the structure and expression levels in a sample of all transcripts that share the same gene locus.
plotMeans('ENSMUSG00000035696.15', bg_filtered, groupvar="group",legend=TRUE)

plotTranscripts(ballgown::geneIDs(bg_filtered)[1729], bg_filtered, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

giurgiu@bioclient10:/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown> ../../Noncoding/softwares/R-3.5.0/bin/R


  x<-subset(df, df$experiment=='103' & df$state == 'B')
    x_entry <- data.frame(FPKM_x = x$FPKM, entry=x$entry, typec=x$typec, replicate=x$replicate)
    y<-subset(df, df$experiment=='103' &  df$state == 'A')
    y_entry <- data.frame(FPKM_y = y$FPKM, entry=y$entry, typec=y$typec, replicate=y$replicate)
    df_aux <- merge(x_entry,y_entry, by.x=c("entry","typec","replicate"), by.y=c("entry","typec","replicate"))
    print(head(df_aux))
    p<- ggplot(df_aux, aes(x=FPKM_x, y=FPKM_y)) +
        geom_point(aes(color = factor(typec))) +
        geom_smooth(method=lm) +
        # scale_y_log10("After",limits=c(1, 1e3)) +
        # scale_x_log10("Before",limits=c(1, 1e3)) +
        ggtitle(paste("M103 FPKM", sep="")) + 
        xlab("Before IP") + ylab("After IP") + 
        theme(legend.position="bottom") + facet_wrap(~factor(replicate),nrow=3)
    
  

------------------------------------------- new ballgown ---------------------------------------------
    
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
    


statisticalTest <- function(sample, size){

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
	write.table(gexpr(bg),file=paste(sample,"_genexpression_fpkm_all.txt", sep=""),sep="\t", quote=F,row.names=T)


############################### filtering


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

####################### PCA Analysis


png(filename=paste("pca_", sample, ".png", sep=""),width=1350, height=850)

########  with FPKM > 1

df_zero = subset(df_genes_expression, df_genes_expression$filterOne == TRUE )
df_zero_transpose = t(df_zero[,1:(2*size)])
data = data.frame(t(df_zero[,1:(2*size)]),"State"=c(rep('After',size),rep('Before',size)))

p1<- autoplot(prcomp(df_zero_transpose), data=data, colour='rownames', shape='State', main=paste("PCA ",sample," (FPKM>1)",sep=""), size=8, xlim=c(-1,1), ylim=c(-1,1))

########  with FPKM > 1 and mean FPKM_a > mean FPKM_b

df_zero = subset(df_genes_expression, df_genes_expression$filterOne == TRUE  & df_genes_expression$filterFPKM_mean == TRUE)
df_zero_transpose = t(df_zero[,1:(2*size)])
data = data.frame(t(df_zero[,1:(2*size)]),"State"=c(rep('After',size),rep('Before',size)))

p2<-autoplot(prcomp(df_zero_transpose), data=data, colour='rownames', shape='State', main=paste("PCA ",sample," (FPKM>1 + mean(PFKM_A)>mean(FPKM_B))",sep=""), size=8, xlim=c(-1,1), ylim=c(-1,1))

########  with FPKM > 0 and mean FPKM_a > mean FPKM_b

df_zero = subset(df_genes_expression, df_genes_expression$filterZero == TRUE  & df_genes_expression$filterFPKM_mean == TRUE)
df_zero_transpose = t(df_zero[,1:(2*size)])
data = data.frame(t(df_zero[,1:(2*size)]),"State"=c(rep('After',size),rep('Before',size)))

p3<-autoplot(prcomp(df_zero_transpose), data=data, colour='rownames', shape='State', main=paste("PCA ",sample," (FPKM>0 + mean(PFKM_A)>mean(FPKM_B))",sep=""), size=8, xlim=c(-1,1), ylim=c(-1,1))

########  with FPKM > 0

df_zero = subset(df_genes_expression, df_genes_expression$filterZero == TRUE )
df_zero_transpose = t(df_zero[,1:(2*size)])
data = data.frame(t(df_zero[,1:(2*size)]),"State"=c(rep('After',size),rep('Before',size)))

p4<-autoplot(prcomp(df_zero_transpose), data=data, colour='rownames', shape='State', main=paste("PCA ",sample," (FPKM>0)",sep=""), size=8, xlim=c(-1,1), ylim=c(-1,1))

########  with mean FPKM_a > mean FPKM_b

df_zero = subset(df_genes_expression, df_genes_expression$filterFPKM_mean == TRUE)
df_zero_transpose = t(df_zero[,1:(2*size)])
data = data.frame(t(df_zero[,1:(2*size)]),"State"=c(rep('After',size),rep('Before',size)))

pca=prcomp(df_zero_transpose)
# pca1
loading_scores_pca1 <- pca$rotation[,1]
genes_scores_pca1 <- abs(loading_scores_pca1)
genes_score_rankes_pca1 <- sort(genes_scores_pca1, decreasing=TRUE)
top_10_pca1 <- names(genes_score_rankes_pca1[1:10])
top_10_pca1 

# pca2
loading_scores_pca2 <- pca$rotation[,2]
genes_scores_pca2 <- abs(loading_scores_pca2)
genes_score_rankes_pca2 <- sort(genes_scores_pca2, decreasing=TRUE)
top_10_pca2 <- names(genes_score_rankes_pca2[1:10])
top_10_pca2

p5<-autoplot(prcomp(df_zero_transpose), data=data, colour='rownames', shape='State', main=paste("PCA ",sample," (mean(PFKM_A)>mean(FPKM_B))",sep=""), size=8, xlim=c(-1,1), ylim=c(-1,1))


########  with geneVariance > 1

df_zero = subset(df_genes_expression, rowVars(df_genes_expression)>1)
df_zero_transpose = t(df_zero[,1:(2*size)])
data = data.frame(t(df_zero[,1:(2*size)]),"State"=c(rep('After',size),rep('Before',size)))

p6<-autoplot(prcomp(df_zero_transpose), data=data, colour='rownames', shape='State', main=paste("PCA ",sample," (var>1)",sep=""), size=8, xlim=c(-1,1), ylim=c(-1,1))

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

dev.off()


}


statisticalTest("M103",5)
statisticalTest("Mlet7",7)



data_stringtie = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/StringTie')
tmp_dir = paste(data_stringtie, sample_IDs, "gene_abundance.tab", sep="/")

transcript_expression = texpr(bg)


######## statistical test with FPKM > 1

	fi=row.names(subset(df_genes_expression, df_genes_expression$filterOne == TRUE & df_genes_expression$filterInfinite == TRUE))

	bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

	# statistical test
	adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)


	adjusted_results = stattest(bg_filtered,feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)



######## statistical test with FPKM > 01and mean FPKM_a > mean FPKM_b

	fi=row.names(subset(df_genes_expression, df_genes_expression$filterOne == TRUE & df_genes_expression$filterFPKM_mean == TRUE))

	bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

	# statistical test
	adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)


	adjusted_results = stattest(bg_filtered,feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)



######## statistical test with FPKM > 0 and mean FPKM_a > mean FPKM_b

	fi=row.names(subset(df_genes_expression, df_genes_expression$filterZero == TRUE & df_genes_expression	$filterFPKM_mean == TRUE))

	bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

	# statistical test
	adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)


	adjusted_results = stattest(bg_filtered,feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)


######## statistical test with FPKM > 0

	fi=row.names(subset(df_genes_expression, df_genes_expression$filterZero == TRUE))

	bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

	# statistical test
	adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)


	# statistical test original
	adjusted_results = stattest(bg_filtered,feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)



######## statistical test mean FPKM_a > mean FPKM_b

	fi=row.names(subset(df_genes_expression, df_genes_expression$filterFPKM_mean == TRUE))

	bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

	# statistical test with matched pairs
	adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
	results_genes = arrange(adjusted_results, qval)
	dim(resuls_genes)
	top14=subset(results_genes, results_genes$qval < 0.05)



	# statistical test original
	adjusted_results = stattest(bg_filtered,feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)


######## statistical test mean FPKM_a > mean FPKM_b + FPKM <= 30000

	fi=row.names(subset(df_genes_expression, df_genes_expression$filterFPKM_mean == TRUE & df_genes_expression$filterInfinite == TRUE ))

	bg_filtered = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

	# statistical test with matched pairs
	adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
	results_genes = arrange(adjusted_results, qval)
	dim(resuls_genes)
	subset(results_genes, results_genes$qval < 0.05)



	# statistical test original
	adjusted_results = stattest(bg_filtered,feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)






######## statistical test variance

	bg_filtered = subset(bg,"rowVars(gexpr(bg)) >1",genomesubset=TRUE)

	# statistical test with matched pairs
	adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
	results_genes = arrange(adjusted_results, qval)
	dim(resuls_genes)
	top13=subset(results_genes, results_genes$qval < 0.05)
	write.table(subset(df_genes_expression,row.names(df_genes_expression) %in% results_genes$id),file="top13_expr_variance.txt",sep="\t",row.names=T,col.names=T, quote=F )


	# statistical test original
	adjusted_results = stattest(bg_filtered,feature='gene',  covariate="group", getFC=TRUE, meas='FPKM')
	results_genes = arrange(adjusted_results, qval)
	subset(results_genes, results_genes$qval < 0.05)

######## filtering isoforms

gene_expression=gexpr(bg)

fpkm = gene_expression
fpkm = log2(fpkm+1)

png(file="M103_geneabundance_new.png")
boxplot(fpkm, col=as.numeric(pData(bg)$group), las=2, ylab='log2(FPKM+1)', main="Distribution of gene abundances across samples (before filtering)")
dev.off

# filter my gene expression
df_genes_expression = data.frame(expression_genes)

# all FPKM > 1
fone=filterOne(expression_genes)
df_genes_expression$filterOne=fone
dim(df_genes_expression[df_genes_expression$filterOne == TRUE,])

# all FPKM > 0
fone=filterZero(expression_genes)
df_genes_expression$filterZero=fone
dim(df_genes_expression[df_genes_expression$filterZero== TRUE,])

######### first filtering

fone = subset(df_genes_expression, df_genes_expression$filterOne == TRUE)
fi = row.names(fone)
bg_filter1 = subset(bg,"ballgown::geneIDs(bg) %in% fi",genomesubset=TRUE)

fpkm = fone[,1:(2*size)]
fpkm = log2(fpkm+1)
png(file="M103_geneabundance_new_fpkm>1.png")
boxplot(fpkm, col=as.numeric(pData(bg_filter1)$group), las=2, ylab='log2(FPKM+1)', main="Distribution of gene abundances across samples (FPKM > 1)")
dev.off()

# count how many LNC have log2(FPKM+1) > 5
dim(subset(fpkm[order(-fpkm$FPKM.M103A1),], fpkm$FPKM.M103A1>5 & gregexpr("^LNC_", row.names(fpkm)) == TRUE))


####################

bg_filter1 = subset(bg,"rowVars(texpr(bg)) >1 & filterHigh(texpr(bg)) == TRUE",genomesubset=TRUE)
gene_ids = unique(ballgown::geneIDs(bg_filter1))

zero = subset(df_genes_expression, df_genes_expression$filterOne== TRUE)

######## second filtering

fi = row.names(subset(zero, row.names(zero) %in% gene_ids))

bg_filtered = subset(bg_filter1,"ballgown::geneIDs(bg_filter1) %in% fi",genomesubset=TRUE)

# statistical test with matched pairs
adjusted_results = stattest(bg_filtered,feature='gene', meas='FPKM', mod0=mod0, mod=mod)
results_genes = arrange(adjusted_results, qval)
head(results_genes)
dim(results_genes)

adjusted_results = stattest(bg_filtered,feature='transcript', meas='FPKM', mod0=mod0, mod=mod)
results_genes = arrange(adjusted_results, qval)
head(results_genes)
dim(results_genes)



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


##############################################################################################################
##############################################################################################################
################################################### CLUSTERING ###############################################
##############################################################################################################
##############################################################################################################



count_t_per_gene <- read.csv2("count_transcripts_per_gene_freq.txt", sep="\t", header=F)
df_transcript2_per_gene <- data.frame(frequency=as.numeric(as.character(count_t_per_gene[,1])), gene=count_t_per_gene[,2], type=count_t_per_gene[,3])



# find genes that need to cluster
df_genes2cluster=subset(df_transcript2_per_gene, df_transcript2_per_gene$frequency > 8 )
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
  
  obj <- collapseTranscripts(gene, bg, meas = "FPKM", method = c("hclust"), k = 4)
  
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

write.table(df,"cluster_tr_8_k_4_M103.txt",sep="\t",row.names=T,quote=F)

# bind FPKMs from other genes that do not 
df_genes=subset(df_transcript2_per_gene, df_transcript2_per_gene$frequency <= 4 )
list_genes=as.character(df_genes$gene)
bg_filter = subset(bg,"ballgown::geneIDs(bg) %in% list_genes",genomesubset=TRUE)
expression_matrix = texpr(bg_filter)
row.names(expression_matrix) <- transcriptNames(bg_filter)
df <- rbind(df, expression_matrix)  
  
# filter matrix with clustered transcripts
df = data.frame(df)

# filtering FPKM > 1
df_filtered = subset(df, filterOne(df) == TRUE)
colnames(df_filtered) <- as.character(pData(bg)$id)
adjusted_results = stattest(gowntable = df_filtered, pData=pData(bg), feature='transcript', meas='FPKM', mod0=mod0, mod=mod)
results_transcript <- arrange(adjusted_results, qval)


#################################################### GET MANUALLY THE GENE EXPRESSION & TRANSCRIPT EXPRESSION
library(plyr)
library(dplyr)

############## read gene expression from file
dir=("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/StringTie")
sample_IDs = list.dirs(path = dir, full.names = FALSE, recursive = FALSE)
sample_IDs <- sample_IDs[ grepl("^M103", sample_IDs) ]
sample_paths = paste(dir,"/", sample_IDs, "/", sample_IDs, "_geneabundance.tab", sep="")

geneIDs_local <- NULL
df <- data.frame(sampleID=character(),
              geneID=character(),
              cov=double(),
              FPKM=double(), 
              TPM=double(), 
              stringsAsFactors=FALSE)

for (i in 1:length(sample_paths)){
  aux <- read.csv2(file=sample_paths[i], sep="\t", header=T)
  df<- rbind(df,data.frame(
                  sampleID=as.character(rep(sample_IDs[i]),dim(aux)[1]),
                  geneID=as.character(aux[,1]),
                  cov=as.numeric(as.character(aux[,7])),
                  FPKM=as.numeric(as.character(aux[,8])),
                  TPM=as.numeric(as.character(aux[,9]))))
  
  print(dim(aux)[1])
  if (is.null(geneIDs_local)){
    geneIDs_local <-  c(as.character(aux[,1]))
  }else{
    geneIDs_new <-  c(geneIDs_local,c(as.character(aux[,1])))
    geneIDs_local <- geneIDs_new
  }
  geneIDs_local <- unique(geneIDs_local)
}

geneIDs_local <- unique(geneIDs_local)

########### GET FPKM
fpkm_matrix <- data.frame(geneID=geneIDs_local)


for (i in 1:length(sample_IDs)){
  replicate<- subset(df, df$sampleID == sample_IDs[i])
  colname<-sample_IDs[i]
  filter_replicate <- data.frame(replicate$FPKM,geneID=replicate$geneID)
  # left join
  fpkm_matrix<-left_join(x=fpkm_matrix, y=filter_replicate , by="geneID",all.x=TRUE,all.y=TRUE)
  colnames(fpkm_matrix)[i+1] <- colname
}

fpkm_matrix<-unique(fpkm_matrix[,1:dim(fpkm_matrix)[2]])

df_aux<-fpkm_matrix[,2:dim(fpkm_matrix)[2]]
row.names(df_aux) <- c(as.character(fpkm_matrix$geneID))
df_filtered = subset(df_aux,filterOne(df_aux) == TRUE)
colnames(df_filtered) <- as.character(pData(bg)$id)
adjusted_results <- stattest(gowntable = df_filtered, pData=pData(bg), feature='gene', meas='FPKM', mod0=mod0, mod=mod)
results_transcript <- arrange(adjusted_results, qval)

########### GET TPM
fpkm_matrix <- data.frame(geneID=geneIDs_local)

df_genes<-data.frame(geneID=geneIDs_local)

for (i in 1:length(sample_IDs)){
    replicate<- subset(df, df$sampleID == sample_IDs[i])
    colname<-sample_IDs[i]
    filter_replicate <- data.frame(replicate$TPM,geneID=replicate$geneID)
    # left join
    fpkm_matrix<-left_join(x=fpkm_matrix, y=filter_replicate , by="geneID")
    colnames(fpkm_matrix)[i+1] <- colname
}

fpkm_matrix<-unique(fpkm_matrix[,1:dim(fpkm_matrix)[2]])
row.names(fpkm_matrix) <- c(as.character(fpkm_matrix$geneID))

df_aux<-fpkm_matrix[,2:dim(fpkm_matrix)[2]]
df_filtered = subset(df_aux,filterOne(df_aux) == TRUE)
colnames(df_filtered) <- as.character(pData(bg)$id)
adjusted_results <- stattest(gowntable = df_filtered, pData=pData(bg), feature='gene', meas='TPM', mod0=mod0, mod=mod)
re <- arrange(adjusted_results, qval)





















