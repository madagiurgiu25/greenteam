
library(ggplot2)
setwd("e:\\masterpraktikum\\atherosclerosis")

library(igraph)
library(d3heatmap)
library(htmlwidgets)
library(UpSetR)

## Read edge list with weights
list1 <- c("mirtrap_mouse_103.txt","mirtrap_human_103.txt","mirtrap_mouse_let7.txt","mirtrap_human_let7.txt")
list2 <- c("plots_intersection_tools/m103_all.png","plots_intersection_tools/h103_all.png","plots_intersection_tools/mlet7_all.png","plots_intersection_tools/hlet7_all.png")
list3 <- c("plots_intersection_tools/m103_lnc.png","plots_intersection_tools/h103_lnc.png","plots_intersection_tools/mlet7_lnc.png","plots_intersection_tools/hlet7_lnc.png")
list4 <- c("plots_intersection_tools/m103_matrix.txt","plots_intersection_tools/h103_matrix.txt","plots_intersection_tools/mlet7_matrix.txt","plots_intersection_tools/hlet7_matrix.txt")

for (i in 1:4){
    data <- read.table(list1[i], header = T, sep = "\t")
    edge_list <- unique(data.frame(V1=data$gene,V2=data$diff_exp_tool))
    
    # edge_list <- subset(edge_list, substr(edge_list$V1,1,3) == 'LNC')
    
    ## Form undirected graph from edge list
    G <- graph.data.frame(edge_list,directed=FALSE)
    ## Get adjacency matrix
    ## Set edge weights to values in the InteractionType column by setting
    ## attr="InteractionType", for an unweighted graph, use attr=NULL
    A<-as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
    
    size_row <- length(unique(edge_list$V1))
    size_col <- length(unique(edge_list$V2))
    
    
    matrix <- A[1:size_row, (size_row + 1):(size_row+size_col)]
    
    m <- data.frame(matrix)
    m[] <- lapply(m, function(x) {
        as.numeric(as.character(x))
    })
    
    png(list2[i],width=800,height=600)
    upset(m,  point.size = 4, line.size = 1, 
          order.by = "freq",
          text.scale = c(2 ,2, 2, 2, 2, 2), sets.x.label = "DiffExp calls per tool", mainbar.y.label = "Databases intersection" )
    dev.off()
    # write.table(m,file=list4[i],row.names = row.names(m), col.names =T, sep="\t",quote=F)
    
    m <- subset(m, substr(row.names(m),1,3) == 'LNC')
    png(list3[i],width=800,height=600)
    upset(m,  point.size = 4, line.size = 1, 
          order.by = "freq",
          text.scale = c(2 ,2, 2, 2, 2, 2), sets.x.label = "DiffExp calls per tool", mainbar.y.label = "Databases intersection" )
    dev.off()
    
    
    # write.table(m,file=list4[i],row.names = row.names(m), col.names =T, sep="\t",quote=F)

}


#### ENA

data <- read.table("ena.txt", header = F, sep = "\t")
edge_list <- unique(data.frame(V1=data$V1,V2=data$V5))

edge_list <- subset(edge_list, substr(edge_list$V1,1,3) == 'ENS')

## Form undirected graph from edge list
G <- graph.data.frame(edge_list,directed=FALSE)
## Get adjacency matrix
## Set edge weights to values in the InteractionType column by setting
## attr="InteractionType", for an unweighted graph, use attr=NULL
A<-as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)

size_row <- length(unique(edge_list$V1))
size_col <- length(unique(edge_list$V2))


matrix <- A[1:size_row, (size_row + 1):(size_row+size_col)]

m <- data.frame(matrix)
m[] <- lapply(m, function(x) {
    as.numeric(as.character(x))
})

png("ena_ens.png",width=800,height=600)
upset(m,  point.size = 4, line.size = 1, 
      order.by = "freq",
      text.scale = c(2 ,2, 2, 2, 2, 2), sets.x.label = "DiffExp calls per Experiment", mainbar.y.label = "LNC expression over the experiments" )
dev.off()
write.table(m,file="ena_matrix_ens.txt",row.names = row.names(m), col.names =T, sep="\t",quote=F)


# write.table(m,file=list4[i],row.names = row.names(m), col.names =T, sep="\t",quote=F)


data <- read.table("plots_intersection_tools/diff_states_il1a_pdgf1.txt", header = F, sep = "\t")

data[order(-V2),]
ggplot(data,aes(x=V1,y=as.numeric(as.character(V2)),group=as.factor(V5),color=as.factor(V5))) + 
    geom_point() + geom_line() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("LNC") + ylab("log2FC") + geom_hline(yintercept=0, linetype="dashed", 
                                              color = "black", size=1) +
    scale_color_discrete(name = "Treatment") + ggtitle("Diff Exp LNC in IL1alpha and PDGF1 (but not in treatment=both)")


setwd("e:\\masterpraktikum\\DiffExp")
datat <- data.frame(read.table("ena.txt.json_protein_coding_coexpression", header = T, sep = "\t"))
datat <- rbind(data,data.frame(read.table("ena.txt.json_miRNA_coexpression", header = T, sep = "\t")))

ggplot(data) + 
geom_point(aes(x=as.numeric(as.character(data$LNC_log2FC)),
               y=as.numeric(as.character(data$ENS_log2FC)))) 
    + facet_wrap( ~ ENS_type) + 
    xlab("LNC_log2FC") + ylab("ENS_log2FC") + 
geom_point(dt = subset(data, data$LNC_tools > 1 & data$ENS_tools > 1),aes(x=as.numeric(as.character(dt$LNC_log2FC)),
                                                                          y=as.numeric(as.character(dt$ENS_log2FC))),colour="red", shape=1, size=2) 

data <- subset(datat,datat$Experiment == 'SRR20546_control_IL1a' & datat$NEIG_type != 'overlap')

with(data, plot(LNC_log2FC, ENS_log2FC, pch=20, main="Coexpression control_IL1a",ylab = "ENS_log2FC",xlab="LNC_log2FC",col="grey",xlim=c(-8,8),ylim=c(-8,8)))
with(subset(data, as.numeric(as.character(data$LNC_count_tools)) > 1 | as.numeric(as.character(data$ENS_count_tools)) > 1), points(LNC_log2FC, ENS_log2FC, pch=20, col="blue"))
with(subset(data, as.numeric(as.character(data$LNC_count_tools)) > 1 & as.numeric(as.character(data$ENS_count_tools)) > 1), points(LNC_log2FC, ENS_log2FC, pch=20, col="red"))
with(data,legend("topleft",c("all", "LNC|ENS (#tools >1)" ,"LNC&ENS (#tools >1)"), fill = c("grey","blue","red")))
with(abline(h=0))
with(abline(v=0))

data <- subset(datat,datat$Experiment == 'SRR20546_control_both' & datat$NEIG_type != 'overlap')

with(data, plot(LNC_log2FC, ENS_log2FC, pch=20, main="Coexpression control_both",ylab = "ENS_log2FC",xlab="LNC_log2FC",col="grey",xlim=c(-8,8),ylim=c(-8,8)))
with(subset(data, as.numeric(as.character(data$LNC_count_tools)) > 1 | as.numeric(as.character(data$ENS_count_tools)) > 1), points(LNC_log2FC, ENS_log2FC, pch=20, col="blue"))
with(subset(data, as.numeric(as.character(data$LNC_count_tools)) > 1 & as.numeric(as.character(data$ENS_count_tools)) > 1), points(LNC_log2FC, ENS_log2FC, pch=20, col="red"))
with(data,legend("topleft",c("all", "LNC|ENS (#tools >1)" ,"LNC&ENS (#tools >1)"), fill = c("grey","blue","red")))
with(abline(h=0))
with(abline(v=0))

data <- subset(datat,datat$Experiment == 'SRR20546_control_PDGF1' & datat$NEIG_type != 'overlap')

with(data, plot(LNC_log2FC, ENS_log2FC, pch=20, main="Coexpression control_PDGF1",ylab = "ENS_log2FC",xlab="LNC_log2FC",col="grey",xlim=c(-8,8),ylim=c(-8,8)))
with(subset(data, as.numeric(as.character(data$LNC_count_tools)) > 1 | as.numeric(as.character(data$ENS_count_tools)) > 1), points(LNC_log2FC, ENS_log2FC, pch=20, col="blue"))
with(subset(data, as.numeric(as.character(data$LNC_count_tools)) > 1 & as.numeric(as.character(data$ENS_count_tools)) > 1), points(LNC_log2FC, ENS_log2FC, pch=20, col="red"))
with(data,legend("topleft",c("all", "LNC|ENS (#tools >1)" ,"LNC&ENS (#tools >1)"), fill = c("grey","blue","red")))
with(abline(h=0))
with(abline(v=0))




###### 
setwd("e:\\masterpraktikum\\atherosclerosis")

data <- read.table("ena_il1a_list.txt", header = F, sep = "\t")
data <- rbind(data,data.frame(read.table("ena_pdgf1_list.txt", header = F, sep = "\t")))
edge_list <- unique(data.frame(V1=data$V1,V2=data$V2))
G <- graph.data.frame(edge_list,directed=FALSE)
A<-as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)

size_row <- length(unique(edge_list$V1))
size_col <- length(unique(edge_list$V2))


matrix <- A[1:size_row, (size_row + 1):(size_row+size_col)]

m <- data.frame(matrix)
m[] <- lapply(m, function(x) {
    as.numeric(as.character(x))
})

heatmap_top100<- d3heatmap(m, scale = "column", colors="RdYlBu", dendrogram="none", Rowv=F)
saveWidget(heatmap_top100, "IL1a_PDGF1.html")


