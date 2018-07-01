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

setwd(
  "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/greenteam/scripts/diffexp/"
)
source("filtering.R", chdir = T)
source("plotsNormalization.R", chdir = T)
source("dataTransformation_library.R", chdir = T)


args = commandArgs(trailingOnly = TRUE)
if (length(args) <= 0) {
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
}

sample = args[1]
size = args[2]
conditions = args[3]
subset_size = args[4]

data_directory = (args[5])
data_bam_directory = (args[6])
working_dir = args[8]

sample = "M103"
size = 5
subset_size = 5
conditions = 2
data_directory = (
  '/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown_Supergenes_exclude20/'
)
data_bam_directory = ('/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/STAR/')



# samples condition
sampleNames_conditions = data.frame(read.csv2(
  file = args[7],
  sep = "\t",
  header = F
))

sampleNames_conditions = data.frame(
  read.csv2(file = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression_mouse/test.txt", sep = "\t", header = F)
)

df_config <- sampleNames_conditions[order(sampleNames_conditions$V1), ]

sampleNames <- as.vector(df_config$V1)
sampleConditions <- as.vector(df_config$V2)

sampleList <- combn(1:size, subset_size, simplify = FALSE)

for (i in 1:length(sampleList)) {
  vector <- unlist(sampleList[i], use.names = F)
  
  # select samples subset
  samples <- lapply(sampleNames, function(x) which(sampleNames == x) %in% vector)
  # name of file for subset
  nameSubset < paste("B", paste(vector, collapse = ''), sep = "")
  
  
}
