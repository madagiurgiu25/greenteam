
library(ggplot2)
library(reshape2)
require("RColorBrewer")
library(gplots)
library(calibrate)


setwd("E:\\masterpraktikum\\mapping")
overlaps <- read.csv("overlaps_formatted.txt",header=T,sep="\t")

lnc_lnc <- subset(overlaps, overlaps$type == 'LNC-LNC')
lnc_ens<- subset(overlaps, overlaps$type == 'LNC-ENS')
ens_ens <- subset(overlaps, overlaps$type == 'ENS-ENS')
lnc_ens_list <- unique(lnc_ens$gene_id_1, lnc_ens$gene_id_2)

lnc_list <- unique(lnc_lnc$gene_id_1, lnc_lnc$gene_id_2)

m <- matrix(0L, nrow = length(lnc_list), ncol = length(lnc_list))
for (i in 1:dim(lnc_lnc)[1]){
    j <- match(lnc_lnc$gene_id_1[i],lnc_list)
    k <- match(lnc_lnc$gene_id_2[i],lnc_list)
    m[j,k] = lnc_lnc$count_overlap[i]
    m[k,j] = lnc_lnc$count_overlap[i]
}

colnames(m) <- lnc_list
row.names(m) <- lnc_list

png(file="heatmap.png")
heatmap.2(m)
dev.off()

########################## M103 expression data analysis
setwd("E:\\masterpraktikum\\mouse")
all<-read.csv2(file="M103_A1235_B1245_expressed_genes.txt",sep="\t",header=T)
significant<-read.csv2(file="M103_A1235_B1245_expressed_genes_significant.txt",sep="\t",header=T)
write.table(row.names(significant),file="M103_A1235_B1245_top_lnc.txt",col.names = F, quote=F,row.names = F)
all$log2FoldChange<-log2(as.numeric(as.character(all$fc)))
all$log10qval <- -log10(as.numeric(as.character(all$qval)))
with(all, plot(log2FoldChange, log10qval, pch=20, main="M103 Diff.exp. genes (A1,2,3,5 - B1,2,4,5)",ylab = "-log10(qval)",col="grey"))
with(all,legend(3.5,2.25,c("lncRNA", "lncRNA, qval < 0.05", "enriched in after IP"), fill = c("black", "blue","red")))
with(subset(all, substring(all$id,1,3) == "LNC" ), points(log2FoldChange, log10qval, pch=20, col="black"))
with(subset(all, as.numeric(as.character(all$qval)) < 0.05 ), points(log2FoldChange, log10qval, pch=20, col="blue"))
with(subset(all, all$id %in% row.names(significant)), points(log2FoldChange, log10qval, pch=20, col="red"))

lnc<-subset(all, substring(all$id,1,3) == "LNC" & as.numeric(as.character(all$qval)) < 0.05  )
write.table(lnc,file="M103_A1235_B1245_top_lnc_values.txt",col.names = T, quote=F,row.names = F,sep="\t")