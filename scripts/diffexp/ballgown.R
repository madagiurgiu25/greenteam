library(ballgown)
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
library(plyr)

# read samples into R
data_directory = ('/Users/Diana/Desktop/greenteam/dataBallgown')
sample_IDs = list.dirs(path = data_directory, full.names = FALSE, recursive = FALSE)
sample_paths = paste(data_directory, sample_IDs, sep="/")
bg = ballgown(samples=sample_paths, meas='all')

# access exon structure e_data
structure(bg)$exon
# access intron structure i_data
structure(bg)$intron
# access transcript structure t_data
structure(bg)$trans

plotTranscripts('NONMMUG034479.2', bg, samples=sample_IDs, meas='FPKM', colorby='transcript')

grep("a", readLines("myfile.dat"), value = TRUE)

pData(bg) = data.frame(id=sampleNames(bg), group=rep(1:0, c(4,4))) # first 4 samples = 1, last 4 samples = 0
phenotype_table = pData(bg)
stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group', getFC=TRUE)
head(stat_results)

# Filter to remove low-abundance genes
bg_filtered = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)

# Identify transcripts that show statistically significant differences between groups.
results_transcripts = stattest(bg_filtered,feature="transcript", covariate="group", libadjust =TRUE, getFC=TRUE, meas="FPKM")

# Identify genes that show statistically significant differences between groups
results_genes = stattest(bg_filtered, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")

# Add gene names to results_transcripts
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filtered),
                                 geneIDs=ballgown::geneIDs(bg_filtered), 
                                 results_transcripts)

# Sort the results from the smallest P value to the largest
results_transcripts = arrange(results_transcripts, pval)
results_genes = arrange(results_genes, pval)

head(results_genes, 10)
head(results_transcripts, 10)


results_transcripts
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
plotMeans('NONMMUG038229.2', bg_filtered, groupvar="group",legend=TRUE)

plotTranscripts(ballgown::geneIDs(bg_filtered)[1729], bg_filtered, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

