
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

edge_list <- subset(edge_list, substr(edge_list$V1,1,3) == 'LNC')

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

png("ena_lnc.png",width=800,height=600)
upset(m,  point.size = 4, line.size = 1, 
      order.by = "freq",
      text.scale = c(2 ,2, 2, 2, 2, 2), sets.x.label = "DiffExp calls per tool", mainbar.y.label = "Databases intersection" )
dev.off()
write.table(m,file="ena_matrix_lnc.txt",row.names = row.names(m), col.names =T, sep="\t",quote=F)

m <- subset(m, substr(row.names(m),1,3) == 'LNC')
png(list3[i],width=800,height=600)
upset(m,  point.size = 4, line.size = 1, 
      order.by = "freq",
      text.scale = c(2 ,2, 2, 2, 2, 2), sets.x.label = "DiffExp calls per tool", mainbar.y.label = "Databases intersection" )
dev.off()


# write.table(m,file=list4[i],row.names = row.names(m), col.names =T, sep="\t",quote=F)


data <- read.table("plots_intersection_tools/diff_states_il1a_pdgf1.txt", header = F, sep = "\t")

data[order(-V2),]
ggplot(data,aes(x=V1,y=as.numeric(as.character(V2)),group=as.factor(V5),color=as.factor(V5))) + 
    geom_point() + geom_line() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("LNC") + ylab("log2FC") + geom_hline(yintercept=0, linetype="dashed", 
                                              color = "black", size=1) +
    scale_color_discrete(name = "Treatment") + ggtitle("Diff Exp LNC in IL1alpha and PDGF1 (but not in treatment=both)")



p<-ggplot(data,aes(x=FID)); 
p+geom_bar(aes(x=factor(FID),y=..count..,fill=STATUS)) 
