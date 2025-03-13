###############################################################################
#'  Created on 13-April-2023
#'  ------------------------------------------------------------
#'  Copyright (c) 2023 Diabetes Epidemiology Group.
#'  All Right Reserved.
#'  ------------------------------------------------------------
#'  Author: anwagbata
#'  Project topic: Scores computation for UKBB
#'  Project info: 
#'     
#'
#'  This is the script for filtering SNPs in the UKBB. The SNPs were
#'  filtered to include SNPs with
#'             (i) minor allele frequency > 0.01
#'            (ii) impute quality info > 0.4
#'            (iii) missingness proportion < 1e-6

library(data.table)

##------------------------------------------------------------------------------
## Loading and using snp-stats-[chr].txt files
##------------------------------------------------------------------------------

setwd("/opt/datastore/shared/ukbb")

## Getting multiple files into R
filelist <- list.files(pattern = "snp-stats-[A-Za-z,0-9]{1,2}\\.txt") 
filelist[1:22] <- gtools::mixedsort(filelist[1:22], decreasing = TRUE)
datalist <- lapply(setNames(, filelist), function(x) fread(x, skip=8L, sep=" ", fill=TRUE))

## Filtering by relevant categories
datalist.filtered <- datalist
datalist.filtered <- lapply(datalist.filtered, function(x) x[x$minor_allele_frequency > 0.01 & 
                                                                          x$impute_info > 0.4 & 
                                                                          abs(x$missing_proportion) < 1E-6 ,])

## Summary of number of SNPs pre and post filter																		
no.row <- lapply(datalist,function(x) nrow(x))
no.row.filtered <- lapply(datalist.filtered,function(x) nrow(x))
pre.filter <- as.data.table(unlist(no.row, recursive = TRUE, use.names = TRUE))
post.filter <- as.data.table(unlist(no.row.filtered, recursive = TRUE, use.names = TRUE))
filter.stat <- cbind(pre.filter, post.filter)

## Rename list objects before saving
new.names <- c("snp-stats-filtered-1","snp-stats-filtered-2","snp-stats-filtered-3",
               "snp-stats-filtered-4","snp-stats-filtered-5","snp-stats-filtered-6",
               "snp-stats-filtered-7","snp-stats-filtered-8","snp-stats-filtered-9",  
               "snp-stats-filtered-10","snp-stats-filtered-11","snp-stats-filtered-12",
               "snp-stats-filtered-13","snp-stats-filtered-14","snp-stats-filtered-15",
               "snp-stats-filtered-16","snp-stats-filtered-17","snp-stats-filtered-18",
               "snp-stats-filtered-19","snp-stats-filtered-20","snp-stats-filtered-21",
               "snp-stats-filtered-22","snp-stats-filtered-X","snp-stats-filtered-XY")
names(datalist.filtered) <- new.names

## Saving files individually
for(i in 1:length(datalist.filtered)){
  write.table(datalist.filtered[[i]], paste0(names(datalist.filtered)[i], ".txt"), 
              row.names=FALSE, quote = FALSE)
}
