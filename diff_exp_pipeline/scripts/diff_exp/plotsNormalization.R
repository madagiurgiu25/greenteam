library(ggplot2)
library(dplyr)
library(stringr)

plotHistogram <- function (df, size, sample, groupconditions, replicates, text){
  
    conditions <- length(unique(groupconditions))
    
    sampledesign <- data.frame(Sample=colnames(df)[1:(conditions * size)],Type=groupconditions,Replicate=replicates)
    
    df_2 <-
        gather(df[, 1:(conditions * size)], "Sample", "FPKM", 1:(conditions * size))
    df_2$FPKM <- as.numeric(as.character(df_2$FPKM))
    df_2 <- merge(df_2, sampledesign)
    # df_2$Type <- substr(df_2$Sample,1,nchar(df_2$Sample)-1)
    # df_2$Replicate <- substr(df_2$Sample, nchar(df_2$Sample), nchar(df_2$Sample))
    
    # p1<-ggplot(df_2, aes(as.numeric(as.character(FPKM)), fill = as.factor(Replicate))) + 
    #   geom_histogram(bins=100, alpha=.5, position="identity") +
    #   xlab("FPKM") + ggtitle(paste("FPKM distribution", sample, text, sep = " ")) +
    #   facet_wrap( ~ Type + Replicate, ncol = conditions) + labs(fill="Replicates")
    
    p1<-ggplot(df_2, aes(as.numeric(as.character(FPKM)))) + 
      geom_histogram(bins=100, alpha=.5, position="identity") +
      xlab("FPKM") + ggtitle(paste("FPKM distribution", sample, text, sep = " ")) +
      facet_wrap( ~ Type + Replicate, ncol = conditions) + labs(fill="Replicates")
    
    
    return(p1)
}

qqplotNormDistr <- function (df, size, sample, conditions){
    
    df_2 <- gather(df[, 1:(conditions * size)], "Sample", "FPKM", 1:(conditions * size))
    df_2$FPKM <- as.numeric(as.character(df_2$FPKM))
    df_2$Type <- substr(df_2$Sample, 1,nchar(df_2$Sample)-1)
    df_2$Replicate <- substr(df_2$Sample, nchar(df_2$Sample), nchar(df_2$Sample))
    
    
    par(mfcol = c(2, size))
    
    for (i in 1:size) {
        data <-
            subset(df_2,
                   df_2$Replicate == i & df_2$Type == paste(sample, "B", sep = ""))
        
        mean_fit <- mean(as.numeric(as.character(data$FPKM)))
        min_fit <- max(as.numeric(as.character(data$FPKM)))
        max_fit <- min(as.numeric(as.character(data$FPKM)))
        sd_fit <- sd(as.numeric(as.character(data$FPKM)))
        
        norm_distr <- rnorm(dim(data$FPKM)[1], mean_fit, sd_fit)    
        
        
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