library(ballgown)
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
    df_trans_aux <- data.frame(trans_data)
    
    return(df_trans_aux)
    
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

transformExp2 <- function(df){
  df[] <- lapply(df, function(x) {
    `^`(2,as.numeric(as.character(x)))
  })
  return(sqrt(df[]))
  
}

quantile_normalization_all <- function(df_trans, size, conditions){
    
    names_col<-colnames(df_trans)
  
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
    colnames(df_norm) <- c("names", "z", names_col[1])
    for (i in 2:(size*conditions)) {
        df_norm <- merge(df_norm, fpkm_list[[i]][1], by.x = "names", by.y = 0)
        names(df_norm)[i + 2] <- names_col[i]
    }
    # for (i in (size+1):(size*conditions)) {
    #     df_norm <- merge(df_norm, fpkm_list[[i]][1], by.x = "names", by.y = 0)
    #     names(df_norm)[i + 2] <- paste(sample, "B", i-size, sep = "")
    # }
    
    row.names(df_norm) <- df_norm$names
    df_norm <- df_norm[,3:(2+conditions*size)]
    
    return(df_norm)
    
}

quantile_normalization_perCondition_new <- function(df_transin,size){
  
  names_col <- colnames(df_transin)
  before_list <- list()
  for (i in 1:size) {
    df <- data.frame(x = as.numeric(as.character(df_transin[,i])))
    row.names(df) <- row.names(df_trans)
    df <- df[order(-df$x), , drop = FALSE]
    before_list[[i]] <- df
  }
  
  for (i in 1:dim(df_transin)[1]) {
    m_b <- 0
    for (j in 1:size) {
      m_b <- m_b + as.numeric(as.character(before_list[[j]][i, ]))
    }
    m_b <- m_b / size
    for (j in 1:size) {
      before_list[[j]][i, ] <- m_b
    }
  }
  
  df_norm1 <- data.frame(z = rep(0, dim(df_transin)[1]))
  row.names(df_norm1) <- row.names(df_transin)
  df_norm1 <- merge(df_norm1, before_list[[1]][1], by = 0)
  colnames(df_norm1) <- c("names", "z", names_col[1])
  for (i in 2:size) {
    df_norm1 <- merge(df_norm1, before_list[[i]][1], by.x = "names", by.y = 0)
    names(df_norm1)[i + 2] <- names_col[i]
  }
  row.names(df_norm1) <- df_norm1$names
  df_norm1 <- df_norm1[,3:dim(df_norm1)[2]]
  
  return(df_norm1)
}

quantile_normalization_perCondition <- function(df_trans,size, conditions){
    
    names_col <- colnames(df_trans)
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
    colnames(df_norm) <- c("names", "z", names_col[1])
    for (i in 2:size) {
        df_norm <- merge(df_norm, after_list[[i]][1], by.x = "names", by.y = 0)
        names(df_norm)[i + 2] <- names_col[i]
    }
    for (i in 1:size) {
        df_norm <- merge(df_norm, before_list[[i]][1], by.x = "names", by.y = 0)
        names(df_norm)[i + 2 + size] <- names_col[i+size]
    }

    row.names(df_norm) <- df_norm$names
    df_norm <- df_norm[,3:(2+conditions*size)]
    
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
createBallgownObj <- function(samplePaths,bamPaths){
  
  bg = ballgown(samples = samplePaths, bamfiles=bamPaths, meas = 'all')
  print(head(ballgown::geneIDs(bg)))
  
  return(bg)
}

callBallgown_geneexpression <- function(bg, df_norm1, subset_size, groupconditions,replicates){
  
  pData(bg) = data.frame(id = sampleNames(bg), group = groupconditions)
  print(pData(bg))
  
  bg_filtered = subset(bg, "ballgown::geneIDs(bg) %in% row.names(df_norm1)", genomesubset = TRUE)
  
  matrix <- as.matrix(mutate_all(df_norm1, function(x) as.numeric(as.character(x))))
  
  colnames(matrix) <- colnames(df_norm1)
  rownames(matrix) <- row.names(df_norm1)
  
  adjusted_results = stattest(
    gowntable = matrix,
    pData=pData(bg_filtered),
    feature = 'gene',
    meas = 'FPKM',
    covariate="group",
    getFC = TRUE,
    libadjust = FALSE
  )
  results_gene <- arrange(adjusted_results, qval)
  
  return(results_gene)
}


callBallgown_2Conditions <- function(bg, df_norm, subset_size){

    pData(bg) = data.frame(id = sampleNames(bg), group = factor(rep(2:1, c(subset_size, subset_size))))
    print(pData(bg))

    bg_filtered = subset(bg,
                         "ballgown::geneIDs(bg) %in% row.names(df_norm)",
                         genomesubset = TRUE)

    matrix <- as.matrix(mutate_all(df_norm, function(x) as.numeric(as.character(x))))

    colnames(matrix) <- colnames(df_norm)
    rownames(matrix) <- row.names(df_norm)

    adjusted_results = stattest(
        gowntable = matrix,
        pData=pData(bg_filtered),
        feature = 'gene',
        meas = 'FPKM',
        covariate="group",
        getFC = TRUE,
        libadjust = FALSE
    )
    results_gene <- arrange(adjusted_results, qval)
    results_filtered <- subset(results_gene, results_gene$qval < 0.05)

    return(results_filtered)
}