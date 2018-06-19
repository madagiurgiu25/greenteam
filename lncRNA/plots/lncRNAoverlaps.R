
library(UpSetR)
# find overlaps between lncRNA


movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )

lnc <- read.csv("E:/masterpraktikum/greenteam/outOverlap.txt",header=T,sep="\t")

lnc_data <- data.frame(lncNamePrimary = lnc$lncNamePrimary,
                       GENCODE = lnc$GENCODE, 
                       NONCODE = lnc$NONCODE, 
                       Lncipedia = lnc$Lncipedia, 
                       lncRNAdb = lnc$lncRNAdb,
                       Exons = lnc$exons_number)

upset(lnc_data,attribute.plots=list(gridrows=60),order.by = "freq", empty.intersections = "on")
      
upset(lnc_data, nsets = 6,  point.size = 4, line.size = 1.5, 
      sets=c("GENCODE","NONCODE","Lncipedia","lncRNAdb"), order.by = "freq",
      queries = list(list(query = intersects, params = list("GENCODE","NONCODE","Lncipedia"), active = T)),
      text.scale = c(1 ,1, 1, 1, 1.5, 1.2), 
      sets.x.label = "lncRNA per database", mainbar.y.label = "Databases intersection")

library(ggplot2)
library(scales)

blank_theme <- theme_minimal()+
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
    )

slices <- c(28469, 1882,56725, 13909)
lbls <- c("lncRNA", "miRNA", "protein-codingRNA", "others")
df <- data.frame(value=slices,type = lbls)
bp<- ggplot(df, aes(x="", y=value, fill=type)) + scale_fill_grey() +  blank_theme +
    theme(axis.text.x=element_blank()) + 
    geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                  label = percent(value/100)), size=5)


pie(slices, labels = lbls, main="GENCODE transcript type - hg38")

slices <- c(17515, 2203,43996, 14121)
lbls <- c("lncRNA", "miRNA", "protein-codingRNA", "others")
pie(slices, labels = lbls, main="GENCODE transcript type - mm10")