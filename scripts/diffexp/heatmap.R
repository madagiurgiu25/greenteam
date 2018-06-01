
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

setwd("E:\\masterpraktikum\\mouse")
data=read.csv("tmp_M103.tab", header=F, sep="\t")
top100_m103=read.csv2("M103_top100genes_info.txt", header=F,sep="\t")
top1000_m103=read.csv2("M103_top1000genes_info.txt", header=F,sep="\t")

m103_tpm <- function(df){
    plots <- list()
    # for (i in 1:4){
    #     print(i)
    #     x<-subset(df, df$experiment=='103' & df$replicate == i & df$state == 'B')
    #     x_entry <- data.frame(TPM_x = x$TPM, entry=x$entry, typec=x$typec)
    #     y<-subset(df, df$experiment=='103' & df$replicate == i & df$state == 'A')
    #     y_entry <- data.frame(TPM_y = y$TPM, entry=y$entry, typec=y$typec)
    #     df_aux <- merge(x_entry,y_entry, by.x=c("entry","typec"), by.y=("entry","typec"))
    #     print(head(df_aux))
    #     p<- ggplot(df_aux, aes(x=TPM_x, y=TPM_y)) +
    #         geom_point(aes(color = factor(typec))) +
    #         geom_smooth(method=lm) +
    #         # scale_y_log10("After",limits=c(1, 1e3)) +
    #         # scale_x_log10("Before",limits=c(1, 1e3)) +
    #         ggtitle(paste("103_",i," TPM", sep="")) + 
    #         xlab("Before IP") + ylab("After IP") + xlim(0,100) + ylim(0,100) +
    #         theme(legend.position="bottom")
    #     plots[[i]] <- p
    #     
    #     
    # }
    
    x<-subset(df, df$experiment=='103' & df$state == 'B')
    x_entry <- data.frame(TPM_x = x$TPM, entry=x$entry, typec=x$typec, replicate=x$replicate)
    y<-subset(df, df$experiment=='103' & df$state == 'A')
    y_entry <- data.frame(TPM_y = y$TPM, entry=y$entry, typec=y$typec, replicate=y$replicate)
    df_aux <- merge(x_entry,y_entry, by.x=c("entry","typec","replicate"), by.y=c("entry","typec","replicate"))
    print(head(df_aux))
    p<- ggplot(df_aux, aes(x=TPM_x, y=TPM_y)) +
        geom_point(aes(color = factor(typec))) +
        geom_smooth(method=lm) +
        ggtitle(paste("M103 TPM", sep="")) + 
        xlab("Before IP") + ylab("After IP") + 
        theme(legend.position="left") + facet_wrap(~factor(replicate),nrow=3)
    
    return(plots)
}

m103_fpkm <- function(df){
    plots <- list()
    # for (i in 1:4){
    #     print(i)
    #     x<-subset(df, df$experiment=='103' & df$replicate == i & df$state == 'B')
    #     x_entry <- data.frame(FPKM_x = x$FPKM, entry=x$entry, typec=x$typec)
    #     y<-subset(df, df$experiment=='103' & df$replicate == i & df$state == 'A')
    #     y_entry <- data.frame(FPKM_y = y$FPKM, entry=y$entry, typec=y$typec)
    #     df_aux <- merge(x_entry,y_entry, by.x=c("entry","typec"), by.y=c("entry","typec"))
    #     print(head(df_aux))
    #     p<- ggplot(df_aux, aes(x=FPKM_x, y=FPKM_y)) +
    #         geom_point(aes(color = factor(typec))) +
    #         geom_smooth(method=lm) +
    #         # scale_y_log10("After",limits=c(1, 1e3)) +
    #         # scale_x_log10("Before",limits=c(1, 1e3)) +
    #         ggtitle(paste("103_",i," FPKM", sep="")) + 
    #         xlab("Before IP") + ylab("After IP") + xlim(0,100) + ylim(0,100) + 
    #         theme(legend.position="bottom") + facet_wrap(~Type,nrow=3)
    #     plots[[i]] <- p
    #     
    # }
    
    x<-subset(df, df$experiment=='103' & df$state == 'B')
    x_entry <- data.frame(FPKM_x = x$FPKM, entry=x$entry, typec=x$typec, replicate=x$replicate)
    y<-subset(df, df$experiment=='103' & df$state == 'A')
    y_entry <- data.frame(FPKM_y = y$FPKM, entry=y$entry, typec=y$typec, replicate=y$replicate)
    df_aux <- merge(x_entry,y_entry, by.x=c("entry","typec","replicate"), by.y=c("entry","typec","replicate"))
    print(head(df_aux))
    p<- ggplot(df_aux, aes(x=FPKM_x, y=FPKM_y)) +
        geom_point(aes(color = factor(typec))) +
        geom_smooth(method=lm) +
        ggtitle(paste("M103 FPKM", sep="")) + 
        xlab("Before IP") + ylab("After IP") + 
        theme(legend.position="left") + facet_wrap(~factor(replicate),nrow=3)
    
    return(p)
}


df_all <- data.frame(state=data$V3, entry = data$V5, sample=data$V1,
                     replicate=data$V4, experiment=data$V2,
                     FPKM=as.numeric(as.character(data$V12)), TPM=as.numeric(as.character(data$V13)),
                     typec=data$V14)

df_103 <- subset(df_all,df_all$experiment=='103')

df <- subset(df_103, entry %in% top100_m103$V2)
df <- df[order(df$entry),]

png(filename = "scatter_fpkm_top100_m103.png", width=600, height = 600)
ml <- marrangeGrob(m103_fpkm(df), nrow=2, ncol=2)
ml
dev.off()

png(filename = "scatter_tpm_top100_m103.png", width=600, height = 600)
ml <- marrangeGrob(m103_tpm(df), nrow=2, ncol=2)
ml
dev.off()

# plot TPM against FPKM

df_aux <- data.frame(FPKM = df$FPKM, entry=df$entry, TPM=df$TPM, state=df$state)
ggplot(df_aux, aes(x=FPKM, y=TPM, color=state)) + geom_point() + scale_y_log10() +
    scale_x_log10() 

df <- subset(df_103, entry %in% top1000_m103$V2)
png(filename = "scatter_fpkm_top1000_m103.png")
ml <- marrangeGrob(m103_fpkm(df), nrow=2, ncol=2)
ml
dev.off()

png(filename = "scatter_tpm_top1000_m103.png")
ml <- marrangeGrob(m103_tpm(df), nrow=2, ncol=2)
ml
dev.off()

df_aux <- data.frame(FPKM = df$FPKM, entry=df$entry, TPM=df$TPM, state=df$state)
ggplot(df_aux, aes(x=FPKM, y=TPM, color=state)) + geom_point() + scale_y_log10() +
    scale_x_log10() 




for (i in 1:4){
    print(i)
    x<-subset(df, df$experiment=='103' & df$replicate == i & df$state == 'B')
    x_entry <- data.frame(FPKM_x = x$FPKM, entry=x$entry)
    y<-subset(df, df$experiment=='103' & df$replicate == i & df$state == 'A')
    y_entry <- data.frame(FPKM_y = y$FPKM, entry=y$entry)
    df_aux <- merge(x_entry,y_entry, by.x="entry", by.y="entry")
    print(head(df_aux))
    p<- ggplot(df_aux, aes(x=FPKM_x, y=FPKM_y)) +
        geom_point(shape=1) +
        geom_smooth(method=lm) +
        # scale_y_log10("After",limits=c(1, 1e5)) +
        # scale_x_log10("Before",limits=c(1, 1e5)) +
        labs(title=paste("103_",i," FPKM", sep=""),xlab="Before",ylab="After")
    plots[[i]] <- p
    
    
}
ml <- marrangeGrob(plots, nrow=2, ncol=2)
ml

# plot TPM against FPKM
x<-subset(df, df$experiment=='103')
df_aux <- data.frame(FPKM = x$FPKM, entry=x$entry, TPM=x$TPM, state=x$state)
ggplot(df_aux, aes(x=FPKM, y=TPM, color=state)) + geom_point() + scale_y_log10() +
    scale_x_log10() 

a<-subset(df, df$experiment== '103' & df$state == 'A' & df$TPM > 100)
b<-subset(df, df$experiment == '103' & df$state == 'B' & df$TPM <10)
extreme<-merge(a,b, by.x=c("entry","replicate"), by.y=c("entry","replicate"))


plots <- list()
count<- 1
for (i in c(1,2,3,4,7,8,9)){
    print(i)
    x<-subset(df, df$experiment=='let' & df$replicate == i & df$state == 'B')
    x_entry <- data.frame(FPKM_x = x$FPKM, entry=x$entry)
    y<-subset(df, df$experiment=='let' & df$replicate == i & df$state == 'A')
    y_entry <- data.frame(FPKM_y = y$FPKM, entry=y$entry)
    df_aux <- merge(x_entry,y_entry, by.x="entry", by.y="entry")
    print(head(df_aux))
    p<- ggplot(df_aux, aes(x=FPKM_x, y=FPKM_y)) +
        geom_point(shape=1) +
        geom_smooth(method=lm) +
        scale_y_log10("After") +
        scale_x_log10("Before") +
        labs(title=paste("let7_",i,sep=""))
    plots[[count]] <- p
    count <- count+1
    
}
ml <- marrangeGrob(plots, nrow=3, ncol=3)
ml

dev.off()

library(dplyr)
# heatmap
library(d3heatmap)
library(htmlwidgets)
df <- subset(df_103, entry %in% top100_m103$V2)

df_new <- data.frame(row.names = top100_m103$V2)
col<- c("M103A1","M103A2","M103A3","M103A4","M103B1","M103B2","M103B3","M103B4")



# M103A1 <- subset(cbind(entry=as.character(df$entry),M103A1=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'A' & df$replicate == 1)
# M103A2 <- subset(cbind(entry=as.character(df$entry),M103A2=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'A' & df$replicate == 2)
# M103A3 <- subset(cbind(entry=as.character(df$entry),M103A3=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'A' & df$replicate == 3)
# M103A4 <- subset(cbind(entry=as.character(df$entry),M103A4=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'A' & df$replicate == 4)
# M103B1 <- subset(cbind(entry=as.character(df$entry),M103B1=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'B' & df$replicate == 1)
# M103B2 <- subset(cbind(entry=as.character(df$entry),M103B2=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'B' & df$replicate == 2)
# M103B3 <- subset(cbind(entry=as.character(df$entry),M103B3=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'B' & df$replicate == 3)
# M103B4 <- subset(cbind(entry=as.character(df$entry),M103B4=as.numeric(as.character(df$FPKM))), df$experiment=='103' & df$state == 'B' & df$replicate == 4)

M103A1 <- subset(cbind(entry=as.character(df$entry),M103A1=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'A' & df$replicate == 1)
M103A2 <- subset(cbind(entry=as.character(df$entry),M103A2=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'A' & df$replicate == 2)
M103A3 <- subset(cbind(entry=as.character(df$entry),M103A3=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'A' & df$replicate == 3)
M103A4 <- subset(cbind(entry=as.character(df$entry),M103A4=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'A' & df$replicate == 4)
M103B1 <- subset(cbind(entry=as.character(df$entry),M103B1=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'B' & df$replicate == 1)
M103B2 <- subset(cbind(entry=as.character(df$entry),M103B2=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'B' & df$replicate == 2)
M103B3 <- subset(cbind(entry=as.character(df$entry),M103B3=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'B' & df$replicate == 3)
M103B4 <- subset(cbind(entry=as.character(df$entry),M103B4=as.numeric(as.character(df$TPM))), df$experiment=='103' & df$state == 'B' & df$replicate == 4)


df_top100 <- data.frame(entry=top100_m103$V2,qval=as.numeric(as.character(top100_m103$V5)))

df_new$entry <- row.names(df_new)
df_new <- merge(df_new,M103A1,by.x="entry",by.y="entry")
df_new <- merge(df_new,M103A2,by.x="entry",by.y="entry")
df_new <- merge(df_new,M103A3,by.x="entry",by.y="entry")
df_new <- merge(df_new,M103A4,by.x="entry",by.y="entry")
df_new <- merge(df_new,M103B1,by.x="entry",by.y="entry")
df_new <- merge(df_new,M103B2,by.x="entry",by.y="entry")
df_new <- merge(df_new,M103B3,by.x="entry",by.y="entry")
df_new <- merge(df_new,M103B4,by.x="entry",by.y="entry")
df_new <- merge(df_new,df_top100,by.x="entry",by.y="entry")
head(df_new)
df_new_sorted<-df_new[order(df_new$qval),]
head(df_new_sorted)
df_new <- df_new_sorted

data <- as.matrix(cbind(as.numeric(as.character(df_new$M103A1)),
                        as.numeric(as.character(df_new$M103A2)),
                        as.numeric(as.character(df_new$M103A3)),
                        as.numeric(as.character(df_new$M103A4)),
                        as.numeric(as.character(df_new$M103B1)),
                        as.numeric(as.character(df_new$M103B2)),
                        as.numeric(as.character(df_new$M103B3)),
                        as.numeric(as.character(df_new$M103B4))))
head(data)
rownames(data) <- df_new$entry
data = na.omit(data)
colnames(data)<-c("M103A1","M103A2","M103A3","M103A4","M103B1","M103B2","M103B3","M103B4")

heatmap_top100<- d3heatmap(data, scale = "column", colors="YlOrRd", dendrogram="none", Rowv=F)
saveWidget(heatmap_top100, "M103_top100.html")
