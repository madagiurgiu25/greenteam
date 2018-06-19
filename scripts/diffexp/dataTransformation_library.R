# library(ballgown)
library(dplyr)

translate_transformlog2 <- function(df, size, conditions){
    
    # convert to numeric
    df[] <- lapply(df, function(x) {
        as.numeric(as.character(x))
    })
    sapply(df, class)
    
    # get min per sample
    minVector <- apply(df[], 2, FUN = min)
    trans_data <- matrix(nrow=dim(df)[1],ncol=dim(df)[2])
    for (j in 1:dim(df)[1]) {
        for (i in 1:(conditions*size)) {
            num <- as.numeric(as.character(df[][j,i]))
            m <- as.numeric(as.character(minVector[[i]]))
            if (num < m){
                print("EROROORSSS")
            }
            trans_data[j,i] = log2(num - m + 1.1)
        }
    }
    row.names(trans_data) <- row.names(df)
    colnames(trans_data) <- colnames(df)
    df_trans <- data.frame(trans_data)
    
    return(df_trans)
    
}


translate_transformlog10 <- function(df, size, conditions){
    
    # convert to numeric
    df[] <- lapply(df, function(x) {
        as.numeric(as.character(x))
    })
    sapply(df, class)
    
    # get min per sample
    minVector <- apply(df[], 2, FUN = min)
    trans_data <- matrix(nrow=dim(df)[1],ncol=dim(df)[2])
    for (j in 1:dim(df)[1]) {
        for (i in 1:(conditions*size)) {
            num <- as.numeric(as.character(df[][j,i]))
            m <- as.numeric(as.character(minVector[[i]]))
            if (num < m){
                print("EROROORSSS")
            }
            trans_data[j,i] = log10(num - m + 1.1)
        }
    }
    row.names(trans_data) <- row.names(df)
    colnames(trans_data) <- colnames(df)
    df_trans <- data.frame(trans_data)
    
    return(df_trans)
    
}

transformSQRT <- function(df, size, conditions){
    # convert to numeric
    df[] <- lapply(df, function(x) {
        as.numeric(as.character(x))
    })
    return(sqrt(df[]))
}

quantile_normalization_all <- function(df_trans, size, conditions){
    
    fpkm_list <- list()
    for (i in 1:(conditions*size)) {
        df <-
            data.frame(x = as.numeric(as.character(df_trans[,i ])))
        row.names(df) <- row.names(df_trans)
        df <- df[order(-df$x), , drop = FALSE]
        fpkm_list[[i]] <- df
    }
    
    for (i in 1:dim(df_trans)[1]) {
        m_b <- 0
        for (j in 1:(conditions*size)) {
            m_b <- m_b + as.numeric(as.character(fpkm_list[[j]][i, ]))
        }
        m_b <- m_b / (conditions*size)
        
        for (j in 1:(conditions*size)) {
            fpkm_list[[j]][i, ] <- m_b
        }
    }
    
    df_norm <- data.frame(z = rep(0, dim(df_trans)[1]))
    row.names(df_norm) <- row.names(df_trans)
    df_norm <- merge(df_norm, fpkm_list[[1]][1], by = 0)
    colnames(df_norm) <- c("names", "z", paste(sample, "A1", sep = ""))
    for (i in 2:size) {
        df_norm <- merge(df_norm, fpkm_list[[i]][1], by.x = "names", by.y = 0)
        names(df_norm)[i + 2] <- paste(sample, "A", i, sep = "")
    }
    for (i in (size+1):(size*conditions)) {
        df_norm <- merge(df_norm, fpkm_list[[i]][1], by.x = "names", by.y = 0)
        names(df_norm)[i + 2] <- paste(sample, "B", i-size, sep = "")
    }
    
    
    row.names(df_norm) <- df_norm$names
    df_norm <- df_norm[,3:(2+conditions*subset_size)]
    
    return(df_norm)
    
}


quantile_normalization_perCondition <- function(df_trans,size, conditions){
    
    before_list <- list()
    after_list <- list()
    for (i in 1:size) {
        df <-
            data.frame(x = as.numeric(as.character(df_trans[,i + size])))
        row.names(df) <- row.names(df_trans)
        df <- df[order(-df$x), , drop = FALSE]
        before_list[[i]] <- df
        
        df <-
            data.frame(x = as.numeric(as.character(df_trans[,i])))
        row.names(df) <- row.names(df_trans)
        df <- df[order(-df$x), , drop = FALSE]
        after_list[[i]] <- df
    }
    
    for (i in 1:dim(df_trans)[1]) {
        m_a <- 0
        m_b <- 0
        for (j in 1:size) {
            m_b <- m_b + as.numeric(as.character(before_list[[j]][i, ]))
            m_a <- m_a + as.numeric(as.character(after_list[[j]][i, ]))
        }
        m_a <- m_a / size
        m_b <- m_b / size
        
        for (j in 1:size) {
            before_list[[j]][i, ] <- m_b
            after_list[[j]][i, ] <- m_a
        }
        
    }
    
    df_norm <- data.frame(z = rep(0, dim(df_trans)[1]))
    row.names(df_norm) <- row.names(df_trans)
    df_norm <- merge(df_norm, after_list[[1]][1], by = 0)
    colnames(df_norm) <- c("names", "z", paste(sample, "A1", sep = ""))
    for (i in 2:size) {
        df_norm <- merge(df_norm, after_list[[i]][1], by.x = "names", by.y = 0)
        names(df_norm)[i + 2] <- paste(sample, "A", i, sep = "")
    }
    for (i in 1:size) {
        df_norm <- merge(df_norm, before_list[[i]][1], by.x = "names", by.y = 0)
        names(df_norm)[i + 2 + size] <- paste(sample, "B", i, sep = "")
    }
    
    row.names(df_norm) <- df_norm$names
    df_norm <- df_norm[,3:(2+conditions*subset_size)]
    
}

variance_stabilization <- function(df, size, conditions){
    
    df[] <- lapply(df, function(x) {
        as.numeric(as.character(x))
    })
    
    for(i in 1:(size*conditions)){
        lambda <- 10/sd(df[][,i])
        df[][,i] <- asinh(df[][,i]/lambda)
    }
    
    return(df[])
}

# callBallgown_2Conditions <- function(bg, df_norm, subset_size){
# 
#     pData(bg) = data.frame(id = sampleNames(bg), group = factor(rep(2:1, c(subset_size, subset_size))))
#     print(pData(bg))
# 
#     bg_filtered = subset(bg,
#                          "ballgown::geneIDs(bg) %in% row.names(df_norm)",
#                          genomesubset = TRUE)
# 
#     matrix <- as.matrix(mutate_all(df_norm, function(x) as.numeric(as.character(x))))
# 
#     colnames(matrix) <- colnames(df_norm)
#     rownames(matrix) <- row.names(df_norm)
# 
#     adjusted_results = stattest(
#         gowntable = matrix,
#         pData=pData(bg_filtered),
#         feature = 'gene',
#         meas = 'FPKM',
#         covariate="group",
#         getFC = TRUE,
#         libadjust = FALSE
#     )
#     results_gene <- arrange(adjusted_results, qval)
#     results_filtered <- subset(results_gene, results_gene$qval < 0.05)
# 
#     return(results_filtered)
# 
# 
# }