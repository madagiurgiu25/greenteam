library(ballgown)
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
library(plyr)
set.seed(1234)

# read samples into R
data_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown')
sample_IDs = list.dirs(path = data_directory, full.names = FALSE, recursive = FALSE)
sample_IDs <- sample_IDs[ grepl("^M103", sample_IDs) ]
sample_paths = paste(data_directory, sample_IDs, sep="/")
bg = ballgown(samples=sample_paths, meas='all')

pData(bg) = data.frame(id=sampleNames(bg), group=rep(1:0, c(4,4))) # first 4 samples = 1, last 4 samples = 0
phenotype_table = pData(bg)
stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group', getFC=TRUE)
#head(stat_results)

# Filter to remove low-abundance genes
bg_filtered = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)

# Identify transcripts that show statistically significant differences between groups.
results_transcripts = stattest(bg_filtered,feature="transcript", covariate="group", getFC=TRUE, meas="FPKM")

# Identify genes that show statistically significant differences between groups
results_genes = stattest(bg_filtered, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")

# Add gene names to results_transcripts
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filtered),
                                 geneIDs=ballgown::geneIDs(bg_filtered), 
                                 results_transcripts)

# Sort the results from the smallest P value to the largest
results_transcripts = arrange(results_transcripts, pval)
results_genes = arrange(results_genes, pval)

head(results_genes, 20)
head(results_transcripts, 20)

write.table(results_genes, file="M103_result_genes_new2.txt", quote=F, sep="\t", row.names = FALSE)
write.table(results_transcripts, file="M103_result_transcripts_new2.txt", quote=F, sep="\t", row.names = FALSE)

tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)

# Show the distribution of gene abundances (measured as FPKM values) across samples, colored by group
fpkm = texpr(bg_filtered, meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm, col=as.numeric(phenotype_table$group), las=2, ylab='log2(FPKM+1)', main="Distribution of gene abundances across samples")

# Plot the structure and expression levels in a sample of all transcripts that share the same gene locus.
plotMeans('ENSMUSG00000035696.15', bg_filtered, groupvar="group",legend=TRUE)

# Identify transcripts and genes with a q value <0.05
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)