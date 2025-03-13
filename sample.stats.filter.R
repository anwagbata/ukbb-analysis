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
#'  This is the script for filtering samples in the UKBB. The samples were
#'  filtered to:
#'      1. Exclude samples with NA values
#'      2. Excluded samples with non-zero values in the last four column:
#'             (i) missing_proportion
#'            (ii) missing_call_proportion
#'            (iii) heterozygous_proportion
#'            (iv) heterozygous_call_proportion

library(data.table)

##------------------------------------------------------------------------------
## Loading and using sample-stats.txt file
##------------------------------------------------------------------------------

sample.stats <- fread("/exports/igmm/eddie/UK-BioBank-Genotype/imputed/v3/sample-stats.txt",
                      header=TRUE, skip = 10L, sep=" ", fill=TRUE)
nrow(sample.stats) 

## Remove samples with NA
sum(is.na(sample.stats))
sample.stats.noNA <- na.omit(sample.stats)
nrow(sample.stats.noNA)

## Remove samples with non-zero in last 4 columns
sample.stats.filtered <- sample.stats.noNA[sample.stats.noNA$missing_proportion ==0 &
                                             sample.stats.noNA$missing_call_proportion ==0 &
                                             sample.stats.noNA$heterozygous_proportion ==0 &
                                             sample.stats.noNA$heterozygous_call_proportion ==0 ,]
nrow(sample.stats.filtered) 

## saving file (${USER})
write.table(sample.stats.filtered,"/exports/eddie/scratch/s1955203/sample-stats-filtered.txt",
            row.names=FALSE, quote = FALSE)
