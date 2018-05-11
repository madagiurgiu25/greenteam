#!/usr/bin/env Rscript

library(limma)

# clear
# cat("\014")  

findDesign <- function (eset){
    
    df_source=data.frame(source=eset$source)
    strings=sort(unique(df_source$source))
    colors=1:length(strings)
    names(colors)=strings
    df_source$source_sort=colors[df_source$source]

    df_description=data.frame(description=eset$description)
    strings=sort(unique(df_description$description))
    colors=1:length(strings)
    names(colors)=strings
    df_description$description_sort=colors[df_description$description]
    
    if (sum(df_source$source_sort) == length(eset$source)){
        return(df_description$description)
    } else if (sum(df_description$description_sort) == length(eset$description)){
        return(df_source$source)
    }else{
        if (max(df_source$source_sort) > max(df_description$description_sort)){
            return(df_source$source)
        }else{
            return(df_description$description)
        }
    }
}

# takes as an input an ExpressionSet (eset), the Design Matrix and the working directory
runDiffExpAnalysis <- function (eset, design){
    fit <- lmFit(eset, design)
    fit <- eBayes(fit)
    topExpressed <- topTable(fit, coef=1, adjust="BH",sort.by = "P", number=Inf)
    return(topExpressed)
    }