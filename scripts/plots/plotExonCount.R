set.seed(1234)
library(ggplot2)
library(reshape2)

plotExonCountHuman <- function(input_file1, input_file2, input_file3, input_file4, separator, title){
  df1 <- read.delim(input_file1, sep = "\t", header=FALSE)
  df2 <- read.delim(input_file2, sep = "\t", header=FALSE)
  df3 <- read.delim(input_file3, sep = "\t", header=FALSE)
  df4 <- read.delim(input_file4, sep = "\t", header=FALSE)
  avg_number_exons1 <- mean(df1$V1)
  print(avg_number_exons1)
  avg_number_exons2 <- mean(df2$V1)
  print(avg_number_exons2)
  avg_number_exons3 <- mean(df3$V1)
  print(avg_number_exons3)
  avg_number_exons4 <- mean(df4$V1)
  print(avg_number_exons4)
  ggplot()+
    geom_density(data = df1, aes(x=df1$V1, y= ..density.., colour="b"), binwidth = 1, stat = "bin", size=1) +
    geom_density(data = df2, aes(x=df2$V1, y= ..density.., colour="r"), binwidth = 1, stat = "bin", size=1) +
    geom_density(data = df3, aes(x=df3$V1, y= ..density.., colour="g"), binwidth = 1, stat = "bin", size=1) +
    geom_density(data = df4, aes(x=df4$V1, y= ..density.., colour="o"), binwidth = 1, stat = "bin", size=1) +
    geom_vline(aes(xintercept=avg_number_exons1), color="black", linetype="dashed") +
    geom_vline(aes(xintercept=avg_number_exons2), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=avg_number_exons3), color="green", linetype="dashed") +
    geom_vline(aes(xintercept=avg_number_exons4), color="orange", linetype="dashed") +
    scale_x_continuous("Number of exons", limits = c(0, 10)) +
    scale_y_continuous("%lncRNA") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(name ="", values=c("r" = "red", "b"="black", "g"="green", "o"="orange"), labels=c("b"="Gencode", "r"="LncRNAdb", "g"="Noncode", "o"="Lncipedia")) +
    scale_fill_manual(values=c("r" = "red", "b"="black", "g"="green", "o"="orange"))+
    theme(legend.position="top", legend.box = "horizontal")
}

plotExonCountMouse <- function(input_file1, input_file2, input_file3, separator, title){
  df1 <- read.delim(input_file1, sep = "\t", header=FALSE)
  df2 <- read.delim(input_file2, sep = "\t", header=FALSE)
  df3 <- read.delim(input_file3, sep = "\t", header=FALSE)
  avg_number_exons1 <- mean(df1$V1)
  print(avg_number_exons1)
  avg_number_exons2 <- mean(df2$V1)
  print(avg_number_exons2)
  avg_number_exons3 <- mean(df3$V1)
  print(avg_number_exons3)
  ggplot()+
    geom_density(data = df1, aes(x=df1$V1, y= ..density.., colour="b"), binwidth = 1, stat = "bin", size=1) +
    geom_density(data = df2, aes(x=df2$V1, y= ..density.., colour="r"), binwidth = 1, stat = "bin", size=1) +
    geom_density(data = df3, aes(x=df3$V1, y= ..density.., colour="g"), binwidth = 1, stat = "bin", size=1) +
    geom_vline(aes(xintercept=avg_number_exons1), color="black", linetype="dashed") +
    geom_vline(aes(xintercept=avg_number_exons2), color="red", linetype="dashed") +
    geom_vline(aes(xintercept=avg_number_exons3), color="green", linetype="dashed") +
    scale_x_continuous("Number of exons", limits = c(0, 10)) +
    scale_y_continuous("%lncRNA") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(name ="", values=c("r" = "red", "b"="black", "g"="green"), labels=c("b"="Gencode", "r"="LncRNAdb", "g"="Noncode")) +
    scale_fill_manual(values=c("r" = "red", "b"="black", "g"="green"))+
    theme(legend.position="top", legend.box = "horizontal")
}

#human
gencode_hg38 = "count_exons_gencode_lnc_hg38_formatted.txt"
lncipedia_hg38 = "count_exons_lncipedia_lnc_hg38_formatted.txt"
lncrnadb_hg38 = "count_exons_lncrnadb_lnc_hg38_formatted.txt"
noncode_hg38 = "count_exons_noncode_lnc_hg38_formatted.txt"

#mouse
gencode_mm10 = "count_exons_gencode_lnc_mm10_formatted.txt"
lncrnadb_mm10 = "count_exons_lncrnadb_lnc_mm10_formatted.txt"
noncode_mm10 = "count_exons_noncode_lnc_mm10_formatted.txt"

plotExonCountHuman(gencode_hg38, lncrnadb_hg38, noncode_hg38, lncipedia_hg38, "\t", "Distribution of number of exons in lncRNA hg38")
plotExonCountMouse(gencode_mm10, lncrnadb_mm10, noncode_mm10, "\t", "Distribution of number of exons in lncRNA mm10")




