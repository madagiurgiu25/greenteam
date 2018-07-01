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
library(here)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("At least one argument must be supplied (input file).n", call. = FALSE)
}

sample = args[1]
size = args[2]
conditions = args[3]

sample = "M103"
size = 5
subset_size = 5
conditions = 2

BallgownPATH = "Ballgown"
TransformationPATH = "Ballgown/Normalization"
data_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown_Supergenes_exclude20/')

setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/greenteam/scripts/diffexp/")
source("filtering.R",chdir=T)
source("plotsNormalization.R",chdir=T)
source("dataTransformation_library.R",chdir=T)

setwd(args[4])

sampleList <- combn(1:size, subset_size, simplify = FALSE)

i = 1

vector <- unlist(sampleList[i], use.names = F)

after <- paste(paste(sample, "A", sep = ""), vector, sep = "")
after_samples <- paste("A", paste(vector, collapse = ''), sep = "")

before <- paste(paste(sample, "B", sep = ""), vector, sep = "")
before_samples <- paste("B", paste(vector, collapse = ''), sep = "")

sample_IDs <- c(after, before)
print(sample_IDs)

# read samples into ballgown object
sample_paths = paste(data_directory, sample_IDs, sep = "/")
print(sample_paths)

bg = ballgown(samples = sample_paths, meas = 'all')
print(head(ballgown::geneIDs(bg)))

# gene expression
expression_genes = gexpr(bg)

# all FPKM > 0
fone = filterZero(expression_genes)
df_genes_expression <- data.frame(expression_genes)
df_genes_expression$filterZero = fone
fmean = filterFPKM_mean(expression_genes)
df_genes_expression$filterFPKM_mean = fmean

# filter data
filtered_data <-
    subset(df_genes_expression, df_genes_expression$filterZero == TRUE)

filtered_data <-
    read.csv2(
        "Ballgown\\filtered_gene_expression_M103.txt",
        sep = "\t",
        row.names = 1,
        header = T
    )

# plot distribution FPKM
plotFPKM <- plotHistogram(filtered_data,subset_size,sample,conditions,"filtered")
plotFPKM

qqplotFPKM <- qqplotNormDistr(filtered_data,subset_size,"FPKM.M103",conditions)
qqplotFPKM


###Translate + Transform log
df_trans <- translate_transformlog2(filtered_data,size,conditions)

plotFPKM <- plotHistogram(df_trans,subset_size,sample,conditions,"transformed")
plotFPKM

qqplotFPKM <- qqplotNormDistr(df_trans,subset_size,"FPKM.M103",conditions)
qqplotFPKM

### Quantile Normalization - overall samples
df_norm <- quantile_normalization_all(df_trans,subset_size, conditions)

plotFPKM <- plotHistogram(df_norm,subset_size,sample,conditions,"transformed")
plotFPKM

qqplotFPKM <- qqplotNormDistr(df_norm,subset_size, sample ,conditions)
qqplotFPKM

########### Quantile Normalization over the 2 conditions
df_norm <- quantile_normalization_perCondition(df_trans,subset_size, conditions)

plotFPKM <- plotHistogram(df_norm,subset_size,sample,conditions,"transformed")
plotFPKM

qqplotFPKM <- qqplotNormDistr(df_norm,subset_size, sample ,conditions)
qqplotFPKM


################### Ballgown call

results <- callBallgown_2Conditions(bg,df_norm, subset_size)
write.table(results,file="M103_translate_logtransform_quantileNorm_diffexp_genes.txt",sep="\t",col.names=T,row.names=T,quote=F)

###################### log10
df_trans <-translate_transformlog10(filtered_data,size,conditions)
df_norm <- quantile_normalization_perCondition(df_trans,subset_size, conditions)

plotFPKM <- plotHistogram(df_norm,subset_size,sample,conditions,"transformed")
plotFPKM

qqplotFPKM <- qqplotNormDistr(df_norm,subset_size, sample ,conditions)
qqplotFPKM

####################################################################################
####################################################################################
#  variance stabilization + Quantile Normalisierung 
df_trans <-variance_stabilization(filtered_data[,1:(size*conditions)],size,conditions)
df_norm <- quantile_normalization_perCondition(df_trans,subset_size, conditions)

plotFPKM <- plotHistogram(df_norm,subset_size,sample,conditions,"transformed")
plotFPKM

qqplotFPKM <- qqplotNormDistr(df_norm,subset_size, sample ,conditions)
qqplotFPKM

####################################################################################
####################################################################################
#  SQRT + Quantile Normalisierung 

###Translate + Transform log
df_trans <- transformSQRT(filtered_data,size,conditions)

plotFPKM <- plotHistogram(df_trans,subset_size,sample,conditions,"transformed")
plotFPKM

qqplotFPKM <- qqplotNormDistr(df_trans,subset_size,"FPKM.M103",conditions)
qqplotFPKM

########### Quantile Normalization over the 2 conditions
df_trans <- translate_transformlog2(filtered_data,size,conditions)
df_norm <- quantile_normalization_perCondition(df_trans,subset_size, conditions)

plotFPKM <- plotHistogram(df_norm,subset_size,sample,conditions,"transformed")
plotFPKM

qqplotFPKM <- qqplotNormDistr(df_norm,subset_size, sample ,conditions)
qqplotFPKM


################### Ballgown call

results <- callBallgown_2Conditions(bg,df_norm, subset_size)



# pData(bg) = data.frame(id = sampleNames(bg), group = factor(rep(2:1, c(
#     subset_size, subset_size
# ))))
# print(pData(bg))
# 
# bg_filtered = subset(bg,
#                      "ballgown::geneIDs(bg) %in% row.names(df_norm)",
#                      genomesubset = TRUE)
# 
# 
matched_pair = factor(rep(1:subset_size, 2))
matched_pair
mod = model.matrix( ~ matched_pair + pData(bg_filtered)$group)
mod0 = model.matrix( ~ pData(bg_filtered)$group)

library(dplyr)
matrix <-
    as.matrix(mutate_all(df_norm, function(x)
        abs(as.numeric(as.character(
            x
        )))))
colnames(matrix) <- colnames(df_norm)
rownames(matrix) <- row.names(df_norm)

adjusted_results = stattest(
    gowntable = matrix,
    pData = pData(bg_filtered),
    feature = 'gene',
    mod0 = mod0,
    mod = mod,
    libadjust = FALSE
)
results_gene <- arrange(adjusted_results, qval)
results_filtered <- subset(results_gene, results_gene$qval < 0.05)
results_filtered
# 
# adjusted_results = stattest(
#     gowntable = matrix,
#     pData=pData(bg_filtered),
#     feature = 'gene',
#     meas = 'FPKM',
#     covariate="group",
#     getFC = TRUE,
#     libadjust = FALSE
# )
# results_gene <- arrange(adjusted_results, qval)
# results_filtered <- subset(results_gene, results_gene$qval < 0.05)
# results_filtered
# 
# genlist <- results_filtered$id
# fpkm_subset <- subset(df_norm, row.names(df_norm) %in% genlist)
# fpkm_subset_original <- subset(filtered_data, row.names(filtered_data) %in% genlist)


# first transformation Quantile Normalization
# CQN (Conditional Quantile Normalization)

library(cqn)
library(scales)
