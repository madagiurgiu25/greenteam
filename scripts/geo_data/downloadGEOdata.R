#!/usr/bin/env Rscript

library(GEOquery)
require("GEOquery")

# setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/GEO")

require(Biobase)

geoID = 'GSE23314'
# geoID = 'GSE28829'
# geoID = 'GSE23304'
geoType = 'GSE' # Dataset, Series or Platform
eset = ''
design = ''

args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  geoID = args[1]
  geoType = args[2]
  pathDir = args[3]
}


# pathDir = "E:\\masterpraktikum\\greenteam\\scripts\\geo_data"
setwd(pathDir)

if (geoType == 'GSE') {
  ############## save the records of interest
  options('download.file.method.GEOquery' = 'libcurl')
  data_matrix <-
    getGEO(geoID,
           GSEMatrix = T,
           AnnotGPL = FALSE,
           destdir =  pathDir)
  
  # save my GSMs (unique ids) that we will use
  for (i in c(1:length(data_matrix))) {
    downloadDetails <-
      pData(phenoData(data_matrix[[i]]))[, c(1, 2, 6, 8, 9, 10, 20, 22, 32)]
    write.table(
      downloadDetails,
      paste(geoID, "_", geoType, "_data.txt", sep = ""),
      col.names = F,
      append = T,
      quote = F,
      sep = "\t"
    )
    
  }
  
  ############## use the data
  data <-
    getGEO(geoID,
           GSEMatrix = F,
           AnnotGPL = FALSE,
           destdir =  pathDir)
  
  ######## Step 1. read meta data
  # get platform for each sample
  gsmplatforms <- do.call('rbind',
                          lapply(GSMList(data), function(x) {
                            Meta(x)$platform_id
                          }))
  # source
  gsmtissue <- do.call('rbind',
                       lapply(GSMList(data), function(x) {
                         Meta(x)$source_name_ch1
                       }))
  # organism
  gsmorganism <- do.call('rbind', lapply(GSMList(data), function(x) {
    Meta(x)$organism_ch1
  }))
  
  # description
  gsmdescription <-
    do.call('rbind', lapply(GSMList(data), function(x) {
      Meta(x)$description
    }))
  
  pdata <- data.frame(
    samples = names(GSMList(data)),
    platform = gsmplatforms,
    organism = gsmorganism,
    source = gsmtissue,
    description = gsmdescription
  )
  
  ######## Step 2. read and format expression matrix
  # make expression matrix
  probesets <- Table(GPLList(data)[[1]])$ID
  
  # row should match the same probeID in the 1..n GSM list
  data.matrix <- do.call('cbind', lapply(GSMList(data), function(x)
  {
    tab <- Table(x)
    mymatch <-
      match(probesets, tab$ID_REF)
    return(tab$VALUE[mymatch])
  }))
  
  data.matrix <-
    apply(data.matrix, 2, function(x) {
      as.numeric(as.character(x))
    })
  data.matrix <- log2(data.matrix)
  
  # create Expression Object
  rownames(data.matrix) <- probesets
  colnames(data.matrix) <- names(GSMList(data))
  
  
  ######## Step 3. build the expression object (matrix + annotation)
  
  # annotations
  pheno <- as(pdata, "AnnotatedDataFrame")
  
  # expression object
  eset <-
    new('ExpressionSet', exprs = data.matrix, phenoData = pheno)
  
}


