

##---------------------------------------------------------------------------

## UKBB genotype directory
/exports/igmm/eddie/UK-BioBank-Genotype/imputed/v3
wc -l /exports/igmm/eddie/UK-BioBank-Genotype//imputed/v3/sample-stats.txt
less /exports/igmm/eddie/UK-BioBank-Genotype//imputed/v3/sample-stats.txt

## Filter sample-stats.txt by excluding individulas with 
library(data.table)
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
nrow(sample.stats.filtered) #475758

## saving file
write.table(sample.stats.filtered,"/exports/eddie/scratch/${USER}/sample-stats-filtered.txt",
            row.names=FALSE, quote = FALSE)


write.table(sample.stats.filtered,"/exports/eddie/scratch/s1955203/snp1.txt", row.names=FALSE, quote = FALSE)


## Copy file to genoscore
rsync -chavzP /exports/eddie/scratch/s1955203/sample-stats-filtered.txt genoscores:/ophead(snp10)t/datastore/shared/ukbb/
  
  ##---------------------------------------------------------------------------  


7402791 = chr1
3388017 = chr10
4931944 = chr3




## On Genoscores Filter snp-stats.txt 


snp1 <- fread("sample-stats-filtered.txt",skip=8L, sep=" ", fill=TRUE)
snp1.filtered <- snp1[snp1$minor_allele_frequency > 0.01 &
                        snp1$impute_info > 0.4 &
                        abs(snp1$missing_proportion) < 1E-6 ,]
nrow(snp1.filtered)

snp1.filtered <- snp1[snp1$minor_allele_frequency > 0.01 ,]

snp1$HW_exact_p_value < 1E-6 &
  
  snp1[HW_exact_p_value < 1E-6 , .N]
snp.hw6 <- snp1[HW_exact_p_value < 1E-6 ,]
snp.hw4 <- snp1[HW_exact_p_value < 1E-4 ,]
snp.hw2 <- snp1[HW_exact_p_value < 1E-2 ,]
snp.maf5 <- snp1[minor_allele_frequency > 0.05 ,]
snp.maf1 <- snp1[minor_allele_frequency > 0.01 ,]

table1 <-  data.table(Var = c(paste0("hw6"),
                              paste0("hw4"),
                              paste0("hw2"),
                              paste0("maf0.05"),
                              paste0("maf0.01")),
                      
                      count = c(paste0(snp1[HW_exact_p_value < 1E-6 , .N]),
                                paste0(snp1[HW_exact_p_value < 1E-4 , .N]),
                                paste0(snp1[HW_exact_p_value < 1E-2 , .N]),
                                paste0(snp1[minor_allele_frequency > 0.05 , .N]),
                                paste0(snp1[minor_allele_frequency > 0.01 , .N])),
                      
                      mean = c(paste0(mean(snp.hw6$HW_exact_p_value)),
                               paste0(mean(snp.hw4$HW_exact_p_value)),
                               paste0(mean(snp.hw2$HW_exact_p_value)),
                               paste0(mean(snp.maf5$minor_allele_frequency)), 
                               paste0(mean(snp.maf1$minor_allele_frequency))),
                      
                      min = c(paste0(min(snp.hw6$HW_exact_p_value)),
                              paste0(min(snp.hw4$HW_exact_p_value)),
                              paste0(min(snp.hw2$HW_exact_p_value)),
                              paste0(min(snp.maf5$minor_allele_frequency)), 
                              paste0(min(snp.maf1$minor_allele_frequency))),
                      
                      max = c(paste0(max(snp.hw6$HW_exact_p_value)),
                              paste0(max(snp.hw4$HW_exact_p_value)),
                              paste0(max(snp.hw2$HW_exact_p_value)),
                              paste0(max(snp.maf5$minor_allele_frequency)), 
                              paste0(max(snp.maf1$minor_allele_frequency))))





##---------------------------------------------------------------------------
## Getting multiple files into R
library(data.table)
setwd("/opt/datastore/shared/ukbb")
filelist <- list.files(pattern = "snp-stats-[A-Za-z,0-9]{1,2}\\.txt") 
filelist[1:22] <- gtools::mixedsort(filelist[1:22], decreasing = TRUE)
datalist <- lapply(setNames(, filelist), function(x) fread(x, skip=8L, sep=" ", fill=TRUE)) 

## Filtering by relevant categories
datalist.filtered <- datalist
datalist.filtered[-23]  <- lapply(datalist.filtered[-23], function(x) x[x$minor_allele_frequency > 0.05 & 
                                                                          x$impute_info > 0.4 & 
                                                                          abs(x$missing_proportion) < 1E-6 ,])

## Chromosome X has different HWE coding categories therefore the filtering above
##didn't work on it, so we filter individually
datalist.filtered[23]  <- lapply(datalist.filtered[23], function(x) x[x$male_female_and_HW_lrt_pvalue < 1E-12 & 
                                                                        x$minor_allele_frequency > 0.05 & 
                                                                        x$impute_info > 0.4 & 
                                                                        abs(x$missing_proportion) < 1E-6 ,])

no.row <- lapply(datalist,function(x) nrow(x))
no.row.filtered <- lapply(datalist.filtered,function(x) nrow(x))

##
pre.filter <- as.data.table(unlist(no.row, recursive = TRUE, use.names = TRUE))
pro.filter <- as.data.table(unlist(no.row.filtered, recursive = TRUE, use.names = TRUE))
cbind(pre.filter, post.filter)


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

##---------------------------------------------------------------------------








ToDo
1. Read QCTools to see if you need to do sample different from SNPs
Tell Andrii x chromosome is 80k
2. Email Michal
3. Qctools, while QCtools is running check epi assessment and materials for UNN conf slide

1. Trans scores for CVD
2. Phenome-wide association for each of our core genes in diabetes so we are 
sure knocking down any gene wont raise BP



filelist

pattern to find something that contains 

2 letter word or number

bash stage.sh sample-stats-filtered.txt





We need to allocate storage on Eddie to store the PLNK files. The exact amount of storage can be inferred only after we convert chromosome 1 to PLINK. Ask Stuart to check what are our options for expanding Eddie storage quota.
Steps to process BGEN files and convert to PLINK files of reasonable size:
  
  Filter BGEN by MAF>0.05 (keep only common SNPs).
Filter by variant imputation quality (qual>0.4).
Filter out all SNPs not in weights table.
Convert BGEN to PLINK using PLINK parallelising on Eddie. Run chromosome 1 to estimate how much storage will be required.
Filter out samples which are not available in our UKBB phenotypes (likely obsolete).




"DONE 1" To filter samples we need to load sample-stats.txt into R data table and 
remove all individuals with non-zeros in last 4 columns and save the updated 
sample-stats.txt in /opt/datastore/shared/ukbb directory on Genoscores as 
sample-stats-filtered.txt. 

"DONE 2" To copy files from Eddie to Genoscores data store you can use the stage.sh script 
in molepi which is run from CLI as: bash stage.sh <filename>. 
(Remember to generate SSH keys on eddie and add them to ~/.ssh/authorized_keys 
  file on Genoscores first).

"DONE 3" To filter SNPs we need to use SNP stats file for each chromosome which are in 
/opt/datastore/shared/ukbb called snp-stats-XX.txt where XX is the chromosome number.
We should read this file in R data.table and filter by:
  
  HW_exact_p_value < 1E-6
minor_allele_frequency > 0.05
impute_info > 0.4
missing_proportion < 3%
save the filtered file to snp-stats-filtered-XX.txt to /opt/datastore/shared/ukbb

"4" Keep only individuals in sample-stats-filtered.txt by running qctools on each 
chromosome BGEN file with -incl-samples flag.

"5" Keep only "good" SNPs in each chromosome by running qctools with flag -incl-rsids.
Filtered BGEN files for each chromosome should be stored on scratch 
/exports/eddie/scratch/<UUI>/, where UUI is your username.

"5" Convert each chromosome to PLINK plink2 --bgen <bgen_file> --make-bed --out <plink_file>, 
where ukb_imp_chrXX is the name of the PLINK file.

On Eddie you could submit job arrays using a script similar to molepi/bgen_subsets.R but 
for your tasks. The analysis is submitted by jobscript.sh. For example in R script you 
could write a qctools command like so:
  cmd <- sprintf("qctool -g %s -snp-stats -osnp %s", bgen.file, stats.file)
system(cmd)
To submit the script for all chromosomes you could do: qsub -t 1-24 jobarray.sh, 
where 1-22 are autosomes and 23 and 24 are sex chromosomes.
Before submitting large jobs such as above I recommend you first try the code in 
interactive session.

The code to run analysis on Eddie is in molepi/ukbb on ukbb branch.
You need to checkout this branch and pull the changes:
  
  
git checkout -b ukbb
git fetch --all
git pull origin/ukbb



$ qctool -g example.bgen -s example.sample -og filtered.bgen -incl-samples samples.txt

This command excludes all samples whose identifier is in the file samples.txt 
(which should contain a whitespace-separated list of identifiers). Samples are 
identified by the first identifier field (often ID_1) in the sample file, or 
if no sample file is specified, by sample identifiers specified in genotype 
data source (e.g. by the header in vcf or bgen formats). The option -incl-samples 
behaves similarly but includes only samples with identifier in the given file. 


$ qctool -g example.bgen -og subsetted.bgen -incl-rsids <filename>
  
  Here the specified file should contain a whitespace-separated list of
rsids that will be excluded from processing


cmd <- sprintf("qctool -g %s -snp-stats -osnp %s", bgen.file, stats.file)
system(cmd)


UKBB genotype data
/exports/igmm/eddie/UK-BioBank-Genotype/imputed/v3

Phenotype in various formats on diabepi server
/opt/shared/project/ukbiobank/data/2023/2023-03-23_V01

information about the variables in dataset - ukb671767.html.

(you can open it with /usr/bin/firefox on the server) and 

information about all the sets of lvariables 
/opt/shared/project/ukbiobank/data/baskets.html



. /etc/profile.d/modules.sh
source /exports/applications/support/set_qlogin_environment.sh
module load igmm/apps/qctool/2.0.8
module load igmm/apps/R/4.1.3
module load R/3.5.3


igmm/apps/bgen/1.1.7                              
module load igmm/apps/bgen/1.2 

-Request a screen session and qlogin
-module load qctools first and see that it works
cp snps filtered filef for chr21 from genoscores to scratch
cp chr21 bgen files to scratch

Go on Bgenix and try excluding 1 SNP
There is option of piping BGEN to QC tool, so if BGEN can exclude SNPs for us the qctool do the sample extraction



## -----------------------------------------------------------------------------
"SNPs filtering"
## -----------------------------------------------------------------------------

library(data.table)
snp21 <- fread("snp-stats-filtered-21.txt")
rsid <- snp21[, c("rsid")]
write.table(rsid, "rsid.txt", row.names=FALSE, quote = FALSE)

## Print file summary for any chromosome
qctool -g ukb_imp_chr21_v3.bgen

################ Using QCTOOLs
screen -S qctooltest
qlogin -l h_vmem=64G -l h_rt=48:00:00
cd /exports/eddie/scratch/s1955203/
module load igmm/apps/qctool/2.0.8
qctool -g ukb_imp_chr21_v3.bgen  -incl-rsids rsid.txt -og ukb_imp_filtered_chr21_v3.bgen
cmd <- sprintf("qctool -g %s -incl-rsids %s -og %s", ukb_imp_chr21_v3.bgen, rsid.txt, ukb_imp_filtered_chr21_v3.bgen )
system(cmd)

################ Using BGENIX
screen -S bgentest
qlogin -l h_vmem=64G -l h_rt=48:00:00
cd /exports/eddie/scratch/s1955203/

## create index file - We already have index file, the bgen.bgi file
bgenix -index -g ukb_imp_chr21_v3.bgen 

## Run
module load igmm/apps/bgen/1.1.7
bgenix -g ukb_imp_chr21_v3.bgen -incl-rsids rsid.txt > outputchr21.bgen

cmd <- sprintf("bgenix -g %s -incl-rsids %s > %s", ukb_imp_chr21_v3.bgen, rsid.txt, outputchr21.bgen)
system(cmd)    

-incl-rsids


## -----------------------------------------------------------------------------
"Sample filtering"
## -----------------------------------------------------------------------------

screen -S plinktest
qlogin -l h_vmem=64G -l h_rt=48:00:00
cd /exports/eddie/scratch/s1955203/
module load igmm/apps/qctool/2.0.8
module load igmm/apps/R/4.1.3
library(data.table)
sample.stat <- fread("sample-stats-filtered.txt")
sample.stats.filtered.id <- sample.stat[, c("sample")]
write.table(sample.stats.filtered.id, "sample-stats-filtered-id.txt", row.names=FALSE, quote = FALSE)


### SAMPLE

qctool -g ukb_imp_chrX_v3.bgen -s ukb23652_imp_chr1_v3_s487371.sample -incl-samples samp_id.txt -og ukb_imp_sampfiltered_chrX_v3.bgen



qctool -g ukb_imp_chr21_v3.bgen -incl-samples testfile.txt -og sample.filtered.bgen 


module load igmm/apps/plink/2.00
plink2 --bgen ukb_imp_chr21_v3.bgen --sample ukb23652_imp_chr1_v3_s487371.sample --make-bed --out chr21

"Problem
Where to find or how to generate .sample file
Test run the sample filtering
Speak to Andrii about all the SNPs having the same number of variants but do wc -l on genoscores before that
View BGEN file of both qctool and bgenix output
Run SNPs in BGENIX and pipe to qctools to filter samples - perform for chromosome 21, will take 2 days to run
for test running, but will need to start looking at GCA in that time
study analysis pipeline in github"


/opt/shared/project/ukbiobank/data/2017/2018-06-12_V03$
  ukb23652_imp_chr1_v3_s487371.sampleZZ

1. Filter Snps using MFI.txt
2. Filter snps using BGENIX
3. Convert to plink
4. Filter samples

qctool 
-g yourdata.bgen 
-s yourdata.sample 
-incl-samples list.txt 
-og yourdata_filtered.bgen

.sample file - Its a type of file in the Oxford format that usually contains 
IDs, missingness, sex, and other covariates (similar to a plinks .fam file).
so in my .sample file do i only need the ids that need to be extracted or all the ids present in the .bgen file?
the order of the ids must be the same as in the genotype file, how can i know the order in the genotype files?
You should have a .sample file together with the BGEN data already. If you don't, it might be worth checking with who provided the data in the first place, since you can technically obtain one from the BGEN file but it's a laborious process.


## -----------------------------------------------------------------------------
"Pipelining QCTOOL and BGENIX"
## -----------------------------------------------------------------------------

As an example, the following command uses bgenix with QCTOOL to compute snp 
summary statistics from a subset of a BGEN file, and view the result on 
the fly using less.

$ bgenix -g file.bgen -range 11:3500000-6500000 | qctool -g - -filetype bgen -snp-stats -osnp - | less -S

###### Piping
bgenix -g ukb_imp_chr21_v3.bgen -incl-rsids rsid.txt | qctool -g ukb_imp_chr21_v3.bgen 
-s yourdata.sample  
-incl-samples sample-stats-filtered-id.txt 
-og sample.filtered.bgen 



bgenix -g ukb_imp_chr21_v3.bgen -incl-rsids rsid.txt > outputchr21.bgen



source("/molepi/shared/process_genotypes/merge.plinkfiles.R")
          
               



cmd <- sprintf("qctool -g %s -snp-stats -osnp %s", bgen.file, stats.file)
system(cmd)

"Quick way to see BGEN file is to convert to VCF using"
bgenix -index -g outputchr21.bgen
bgenix -g outputchr21.bgen -vcf




Convert to Plink
Note that PLINK 2 collapses the raw probabilities stored in .gen/.bgen files 
down to dosages; you cannot use PLINK 2 to losslessly convert between 
e.g. BGEN sub-formats. (But if the next program in your 

plink2 --bgen data.chr22.bgen --sample [filename] --make-bed --out ex22

plink2 \
--bgen $imputeddata/imp_data_chr${1}_maf1_info4.bgen ref-first \
--sample $imputeddata/data_chr1.sample \
--extract $resources/w_hm3.snplist \
--make-just-bim \
--out $projectdata/ukb_imp_hm3_chr${1}

## Removing SNPs with filtering categories as NA
snp21 <- fread("/opt/datastore/shared/ukbb/snp-stats-21.txt", skip=8L, sep=" ", fill=TRUE)

## Filtering by relevant categories
snp21.filtered <- snp21[snp21$HW_exact_p_value < 1E-6 &
                          snp21$minor_allele_frequency > 0.05 &
                          snp21$impute_info > 0.4 &
                          abs(snp21$missing_proportion) < 1E-6 ,]
nrow(snp21.filtered)
datalist.nona <- lapply(datalist, function(x) x[!is.na(x$HW_exact_p_value) & 
                                                  !is.na(x$minor_allele_frequency) &
                                                  !is.na(x$impute_info) &
                                                  !is.na(x$missing_proportion) , ]
                        
                        
                        ## Removing SNPs with relevant filtering categories as NA
                        snp21 <- snp21[!is.na(snp21$HW_exact_p_value) &
                                         !is.na(snp21$minor_allele_frequency) &
                                         !is.na(snp21$impute_info) &
                                         !is.na(snp21$missing_proportion) , ]
                        
                        
                        
                        - bgen(v1.2): genetic dosage in binary format include rsid or chromosome_position_ref_alt
                        - bgi: bgen index similar to tbi/cbi tabix indexes
                        - mfi: variant informations: rsid, position (bp), reference allele, alternate allele, maf, info(imputation quality)
                        
                        
                        
Passing arguments to an R script from command lines
https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
https://stackoverflow.com/questions/2151212/how-can-i-read-command-line-parameters-from-an-r-script






module load igmm/apps/plink/2.00
plink2 --bgen ukb_imp_chr21_v3.bgen --sample ukb23652_imp_chr1_v3_s487371.sample --make-bed --out chr21

bgen.file <- file.path(bgen.dir, sprintf("ukb_imp_chr%s_v3.bgen", chr))
stats.file <- file.path(out.dir, sprintf("ukb_mfi_filtered_chr%s_v3.txt", chr))
output.file <- file.path(out.dir, sprintf("ukb_imp_filtered_chr%s_v3.bgen", chr))

cmd <- sprintf("bgenix -g %s -include-rsids %s > %s", bgen.file, stats.file, out.file)
system(cmd)




                        