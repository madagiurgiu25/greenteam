#!/usr/bin/env Rscript

library(GEOquery)
setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/GEO")
require(Biobase)

geoID = 'GSE23314'
geoType = 'GSE' # Dataset, Series or Platform
eset = ''
design = ''

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
	stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else {
	geoID = args[1]
	geoType = args[2]
} 

if (geoType == 'GSE'){

	############## save the records of interest
	data_matrix <- getGEO(geoID,GSEMatrix=T)
	
	# save my GSMs (unique ids) that we will use 
	for(i in c(1:length(data))){
		downloadDetails<-pData(phenoData(data[[i]]))[,c(1,2,6,8,9,10,20,22,32)])
		write.table(downloadDetails, paste(geoID,"_",geoType,"_data.txt", sep = "\t", col.names = F, append = T, quote = F)
		
	}

	############## use the data
	data <- getGEO(geoID,GSEMatrix=F)
	# get platform for each sample
	gsmplatforms <- lapply(GSMList(data),function(x) {Meta(x)$platform})
	gsmtissue <- lapply(GSMList(data),function(x) {Meta(x)$source_name_ch1})
	gsmorganism <- lapply(GSMList(data),function(x) {Meta(x)$organism})

	probesets <- Table(GPLList(data)[[1]])$ID

	# row should match the same probeID in the 1..n GSM list
	data.matrix <- do.call('cbind', lapply(GSMList(data),function(x)
                                       {tab <- Table(x)
                                        mymatch <- match(probesets,tab$ID_REF)
                                        return(tab$VALUE[mymatch])
                                      }))

	data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
	data.matrix <- log2(data.matrix)
	
	# create Expression Object
	rownames(data.matrix) <- probesets
	colnames(data.matrix) <- names(GSMList(data))
	pdata <- data.frame(samples=names(GSMList(data)),organism=gsmorganism,tissue=gsmtissue)
	
	
	pdata <- do.call('rbind',lapply(GSMList(data),function(x)
                                       {org <- Meta(x)$organism
										tissue <- Meta(x)$source_name_ch1
                                        sample <- Meta(x)$geo_accession
                                        return(c(org,tissue,sample))
                                      }))
	pdataFrame <- data.frame(samples = pdata[,3],
							tissue = pdata[,2],
							organism = pdata[,1])
	
	rownames(pdata) <- names(GSMList(data))
	
	pheno <- as(pdata,"AnnotatedDataFrame")

	# expression object
	eset <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)

}

return eset
