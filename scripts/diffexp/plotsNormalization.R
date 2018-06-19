library(ggplot2)
library(dplyr)
library(stringr)

plotHistogram <- function (df, size, sample, conditions,text){
    df_2 <-
        gather(df[, 1:(conditions * size)], "Sample", "FPKM", 1:(conditions * size))
    df_2$FPKM <- as.numeric(as.character(df_2$FPKM))
    df_2$Type <- substr(df_2$Sample,1,nchar(df_2$Sample)-1)
    df_2$Replicate <- substr(df_2$Sample, nchar(df_2$Sample), nchar(df_2$Sample))
    
    p1<- ggplot(df_2, aes(as.numeric(as.character(FPKM)), fill = Replicate)) +
        geom_histogram(bins = 100) +
        xlab("FPKM") + ggtitle(paste("FPKM distribution", sample, text, sep = " ")) +
        facet_wrap( ~ Type, ncol = conditions)
    return(p1)
}

qqplotNormDistr <- function (df, size, sample, conditions){
    
    df_2 <- gather(df[, 1:(conditions * size)], "Sample", "FPKM", 1:(conditions * size))
    df_2$FPKM <- as.numeric(as.character(df_2$FPKM))
    df_2$Type <- substr(df_2$Sample, 1,nchar(df_2$Sample)-1)
    df_2$Replicate <- substr(df_2$Sample, nchar(df_2$Sample), nchar(df_2$Sample))
    
    
    par(mfcol = c(2, size))
    mean_fit <- as.numeric(as.character(mean(df_2$FPKM)))
    min_fit <- as.numeric(as.character(min(df_2$FPKM)))
    max_fit <- as.numeric(as.character(max(df_2$FPKM)))
    
    norm_distr <- rnorm(dim(df)[1], mean_fit, 1)
    print(head(norm_distr))
    
    for (i in 1:size) {
        data <-
            subset(df_2,
                   df_2$Replicate == i & df_2$Type == paste(sample, "B", sep = ""))
        qqplot(
            quantile(norm_distr, probs = seq(0.01, 0.99, 0.01)),
                y = as.numeric(as.character(data$FPKM)),
            ylim = c(min_fit, max_fit),
            xlab = "Theoretical Quantiles - normal distribution",
            ylab = "Sample Quantiles (Before)",
            main = paste("Normal Q-Q Plot ", sample, "_", toString(i), sep = "")
        ); qqline(y = as.numeric(as.character(data$FPKM)), col = 2,lwd=2,lty=2)
        
        data <-
            subset(df_2,
                   df_2$Replicate == i & df_2$Type == paste(sample, "A", sep = ""))
       qqplot(
            quantile(norm_distr, probs = seq(0.01, 0.99, 0.01)),
            y = as.numeric(as.character(data$FPKM)),
            ylim = c(min_fit, max_fit),
            xlab = "Theoretical Quantiles - normal distribution",
            ylab = "Sample Quantiles (After)",
            main = paste("Normal Q-Q Plot ", sample, "_", toString(i), sep = "")
        ); qqline(y = as.numeric(as.character(data$FPKM)), col = 2,lwd=2,lty=2)
        
        
    }
}