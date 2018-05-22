# plot length distribution and overlay between 2 datasets, column 'length'
plotLengthDistribution <- function(input_file1, input_file2, separator, title){
  df1 <- read.delim(input_file1, sep = "\t")
  df2 <- read.delim(input_file2, sep = "\t")
  avg_length1 <- mean(df1$length)
  print(avg_length1)
  avg_length2 <- mean(df2$length)
  print(avg_length2)
  ggplot()+
    geom_histogram(data = df1, aes(x=df1$length, fill="b", colour="b"), binwidth = 5, fill="white", alpha=0.5, stat = "bin") +
    geom_histogram(data = df2, aes(x=df2$length, fill="r", colour="r"), binwidth = 5, fill="white", alpha=0.5, stat = "bin") +
    geom_vline(aes(xintercept=avg_length1), color="black", linetype="dashed", size = 1) +
    geom_vline(aes(xintercept=avg_length2), color="red", linetype="dashed", size = 1) +
    scale_x_continuous("miRNA length") +
    scale_y_continuous("Frequency") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(name ="", values=c("r" = "red", "b"="black"), labels=c("b"="hchg38", "r"="hcmm10")) +
    scale_fill_manual(values=c("r" = "red", "b"="black"))+
    theme(legend.position="top", legend.box = "horizontal")
}

set.seed(1234)
library(ggplot2)
library(reshape2)
intersect_hg38_lengths = "intersect_hg38_lengths.txt"
intersect_mm10_lengths = "intersect_mm10_lengths.txt"
plotLengthDistribution(intersect_hg38_lengths, intersect_mm10_lengths, "\t", "miRNA length distribution in hchg38 and hcmm10")

mirbase_mi_hg38_lengths = "mirbase_mi_hg38_lengths.txt"
gencode_mi_hg38_lengths = "gencode_mi_hg38_lengths.txt"
plotLengthDistribution(mirbase_mi_hg38_lengths, gencode_mi_hg38_lengths, "\t", "miRNA length distribution in hg38")

mirbase_mi_mm10_lengths = "mirbase_mi_mm10_lengths.txt"
gencode_mi_mm10_lengths = "gencode_mi_mm10_lengths.txt"
plotLengthDistribution(mirbase_mi_mm10_lengths, gencode_mi_mm10_lengths, "\t", "miRNA length distribution in mm10")