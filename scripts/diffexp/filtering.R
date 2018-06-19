
filterOne <- function(expr_data, na.rm = TRUE){
    size <- dim(expr_data)[2]
    results <- c(rep(TRUE,dim(expr_data)[1]))
    for (i in 1:dim(expr_data)[1]){
        for(j in 1:size)
            if(expr_data[i,j] < 1 ){
                results[i] <- FALSE
            }
    }
    return(results)
}

filterInfinite <- function(expr_data, na.rm = TRUE){
    size <- dim(expr_data)[2]
    results <- c(rep(TRUE,dim(expr_data)[1]))
    for (i in 1:dim(expr_data)[1]){
        for(j in 1:size)
            if(expr_data[i,j] > 10000 ){
                results[i] <- FALSE
            }
    }
    return(results)
}


filterHigh <- function(expr_data, na.rm = TRUE){
    size <- dim(expr_data)[2]
    results <- c(rep(TRUE,dim(expr_data)[1]))
    for (i in 1:dim(expr_data)[1]){
        for(j in 1:size)
            if(expr_data[i,j] > 500 ){
                results[i] <- FALSE
            }
    }
    return(results)
}


filterZero <- function(expr_data,  na.rm = TRUE){
    size <- dim(expr_data)[2]
    results <- c(rep(TRUE,dim(expr_data)[1]))
    for (i in 1:dim(expr_data)[1]){
        for(j in 1:size)
            if(expr_data[i,j] ==0 ){
                results[i] <- FALSE
            }
    }
    return(results)
}

filterZero2 <- function(expr_data, na.rm = TRUE){
    size <- dim(expr_data)[2]
    results <- c(rep(TRUE,dim(expr_data)[1]))
    for (i in 1:dim(expr_data)[1]){
        count <- 0
        for(j in 1:(size/2))
            if(expr_data[i,j] !=0 ){
                count <- count +1
            }
        if (count< (size/2 -1)){
            results[i] <- FALSE
        }else{
            count <- 0
            for(j in (size/2+1):size)
                if(expr_data[i,j] !=0 ){
                    count <- count +1
                }
            if (count< (size/2 -1)){
                results[i] <- FALSE
            }
        }
    }
    return(results)
}

filterFPKM <- function(expr_data, na.rm = TRUE){
    size <- dim(expr_data)[2]
    results <- c(rep(TRUE,dim(expr_data)[1]))
    for (i in 1:dim(expr_data)[1]){
        for(j in 1:(size/2))
            if(expr_data[i,j] < expr_data[i,size/2 - j + 2] ){
                results[i] <- FALSE
            }
    }
    return(results)
}

filterFPKM_mean <- function(expr_data, na.rm = TRUE){
    size <- dim(expr_data)[2]
    results <- c(rep(TRUE,dim(expr_data)[1]))
    for (i in 1:dim(expr_data)[1]){
        half = size/2
        after = sum(expr_data[i,1:half])/half
        
        before = sum(expr_data[i,(half+1):size])/half
        if (before > after){
            results[i] <- FALSE
        }
        
    }
    return(results)
}
