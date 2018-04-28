# update bioconductor to 3.6
source("https://bioconductor.org/biocLite.R")
biocLite()

# get GEOquery
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
biocLite(c("Biobase"))
biocLite("limma")

library(GEOquery)
library(DESeq)
library(Biobase)
library(limma)

setwd("/home/g/giurgiu/Dokumente/master_sem2/Masterpraktikum/geo_data")
file <- "geo_accessions.txt"
destDir <- "/home/g/giurgiu/Dokumente/master_sem2/Masterpraktikum/geo_data"
geoObj <- getGEO(GEO = 'GDS4527', filename = NULL, destdir = destDir,
                 GSElimits =NULL, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE,
                 parseCharacteristics = TRUE)
# expression set
eset <- GDS2eSet(geoObj,do.log2=TRUE)

gpl_name <- Meta(geoObj)$platform
gpl <- getGEO(GEO = gpl_name, filename = NULL, destdir = destDir,
              GSElimits =NULL, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE,
              parseCharacteristics = TRUE)

## ------------------------------------------------------------------------
# expression set
MA <- GDS2MA(geoObj,GPL=gpl)
class(MA)

design <- model.matrix(~ 0+factor(c(rep.int(1,10),rep.int(2,10))))
colnames(design) <- unique(MA$targets[[2]])
fit <- lmFit(eset, design)
fit <- eBayes(fit)
topTable(fit, coef=1, adjust="BH")
write.table(fit,"GDS4527_diff_exp.txt")

arrayw <- arrayWeights(MAlms)
barplot(arrayw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)

# fit <- lmFit(MA)
# fit <- eBayes(fit)
# boxplot(fit$t~MA$genes$Status, at=1:5-0.2, col=5, boxwex=0.4, xlab="control type",
#           ylab="moderated t-statistics", pch=".", ylim=c(-70, 70), medlwd=1)
# boxplot(fitw$t~MA$genes$Status, at=1:5+0.2, col=6, boxwex=0.4,
#           add=TRUE, yaxt="n", xaxt="n", medlwd=1, pch=".")
# abline(h=0, col="black", lty=2, lwd=1)
# legend(0.5, 70, legend=c("Equal weights", "Array weights"), fill=c(5,6), cex=0.8)

