#!/usr/bin/env Rscript

library(GEOquery)
require(GEOquery)
require(Biobase)

analyzeGEO <- function(geoType, geoID, pathDir) {
    options('download.file.method.GEOquery' = 'libcurl')
    
    if (geoType == 'GSE') {
        data_info <-
            getGEO(
                geoID,
                GSEMatrix = T,
                AnnotGPL = T,
                destdir =  pathDir
            )
        
        for (i in c(1:length(data_info))) {
            # write GSM and pheno data
            ############## save the records of interest
            downloadDetails <-
                pData(phenoData(data_info[[i]]))[, c(1, 2, 6, 8, 9, 10, 20, 22, 32)]
            write.table(
                downloadDetails,
                paste(geoID, "_", geoType, "_pData_", i, ".txt", sep = ""),
                col.names = F,
                append = T,
                quote = F,
                sep = "\t"
            )
            
            # get platform for each sample
            gsmplatforms <- do.call('rbind',as.list(as.vector(data_info[[i]]@phenoData@data$platform_id)))
            
            # source
            gsmtissue <- do.call('rbind',as.list(as.vector(data_info[[i]]@phenoData@data$source_name_ch1)))
            
            # organism
            gsmorganism <- do.call('rbind',as.list(as.vector(data_info[[i]]@phenoData@data$organism_ch1)))
            
            # description
            gsmdescription <- do.call('rbind',as.list(as.vector(data_info[[i]]@phenoData@data$description)))
            
            pdata <- data.frame(
                samples = sampleNames(data_info[[i]]),
                platform = gsmplatforms,
                organism = gsmorganism,
                source = gsmtissue,
                description = gsmdescription
            )
            row.names(pdata) <- sampleNames(data_info[[i]])
            
            ######## Step 2. read and format expression matrix
            # make expression matrix
            probesets <- data_info[[i]]@featureData@data$ID
            
            data.matrix <- exprs(data_info[[i]])
            
            data.matrix <-
                apply(data.matrix, 2, function(x) {
                    as.numeric(as.character(x))
                })
            
            data.matrix <- log2(data.matrix)
     
            ######## Step 3. build the expression object (matrix + annotation)
            
            # annotations
            pheno <- as(pdata, "AnnotatedDataFrame")
            rownames(data.matrix) <- probesets
            colnames(data.matrix) <- sampleNames(data_info[[i]])
            
            # expression object
            eset <- 
                new('ExpressionSet',
                    exprs = data.matrix,
                    phenoData = pheno)
            
            ############## use the data
            # data <-
            #     getGEO(
            #         geoID,
            #         GSEMatrix = F,
            #         AnnotGPL = T,
            #         destdir =  pathDir
            #     )
            # Print the Details for each experiment (sometimes a series contains multiple experiments)
            # show(pData(data[[i]]))
            
            
            ######## Step 1. read meta data
            # get platform for each sample
            # gsmplatforms <- do.call('rbind',
            #                         lapply(GSMList(data), function(x) {
            #                             Meta(x)$platform_id
            #                         }))
            # # source
            # gsmtissue <- do.call('rbind',
            #                      lapply(GSMList(data), function(x) {
            #                          Meta(x)$source_name_ch1
            #                      }))
            # # organism
            # gsmorganism <-
            #     do.call('rbind', lapply(GSMList(data), function(x) {
            #         Meta(x)$organism_ch1
            #     }))
            # 
            # # description
            # gsmdescription <-
            #     do.call('rbind', lapply(GSMList(data), function(x) {
            #         Meta(x)$description
            #     }))
            # 
            # pdata <- data.frame(
            #     samples = names(GSMList(data)),
            #     platform = gsmplatforms,
            #     organism = gsmorganism,
            #     source = gsmtissue,
            #     description = gsmdescription
            # )
            # 
            ######## Step 2. read and format expression matrix
            # make expression matrix
            # probesets <- Table(GPLList(data)[[1]])$ID
            
            # row should match the same probeID in the 1..n GSM list
            # data.matrix <-
            #     do.call('cbind', lapply(GSMList(data), function(x)
            #     {
            #         tab <- Table(x)
            #         mymatch <-
            #             match(probesets, tab$ID_REF)
            #         return(tab$VALUE[mymatch])
            #     }))
            # 
            # data.matrix <-
            #     apply(data.matrix, 2, function(x) {
            #         as.numeric(as.character(x))
            #     })
            # data.matrix <- log2(data.matrix)
            
            # create Expression Object
            # rownames(data.matrix) <- probesets
            # colnames(data.matrix) <- names(GSMList(data))
            
            
            ######## Step 3. build the expression object (matrix + annotation)
            
            # annotations
            # pheno <- as(pdata, "AnnotatedDataFrame")
            # 
            # # expression object
            # eset <-
            #     new('ExpressionSet',
            #         exprs = data.matrix,
            #         phenoData = pheno)
            
            
            
            ####### Step 4. run differential expression analysis
            runDiffExpAnalysis(
                eset = eset,
                design = c(1, 1, 1, 1),
                wdir = pathDir
            )
        }
    }
}



# geoID = 'GSE23314'
geoID = 'GSE28829'
geoID = 'GSE4840'
# geoID = 'GSE22543'
# geoID = 'GSE23304'
geoType = 'GSE' # Dataset, Series or Platform
eset = ''
design = ''

args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
    stop("At least one argument must be supplied (input file).\n",
         call. = FALSE)
} else {
    geoID = args[1]
    geoType = args[2]
    pathDir = args[3]
}


pathDir = "E:\\masterpraktikum\\greenteam\\scripts\\geo_data"
setwd(pathDir)
source("E:\\masterpraktikum\\greenteam\\scripts\\geo_data\\diffExpressionAnalysis.R")

# read GEO entry and make the diff expression analysis
analyzeGEO(geoType = geoType,
           geoID = geoID,
           pathDir = pathDir)

# get the features names (in our case the expressed transcripts)
ids <- featureNames(eset)

# expression matrix rows = features, columns = samples
matrix <- exprs(eset)

design <- model.matrix(~ 0+factor(cbind(df_source$source_sort)))
colnames(design) <- cbind("krank","normal")
fit <- lmFit(eset, design)
fit <- eBayes(fit)
topTable(fit, coef=1, adjust="BH")
