rm(list=ls())
library(tidyverse)
cellinfo <- data.table::fread("../E-MTAB-5061.sdrf.txt") %>%
    select(id=`Comment[ENA_RUN]`, cell=`Source Name`, 
           type=`Characteristics[cell type]`, subj=`Characteristics[individual]`) %>%
    filter(type!="not applicable")
subjvec <- unique(cellinfo$subj)

temp <- as.integer(commandArgs(trailingOnly = TRUE))
donor <- subjvec[temp]

cellinfo.subj <- cellinfo %>% filter(subj==donor)
picard.code <- paste("picard MergeSamFiles",
                     paste(paste0("I=../bam_trimmed_filtered/",cellinfo.subj$id,"Aligned.sortedByCoord.markdup.SplitNCigar.readgrps.BQSR.bam"),collapse=" "),
                     paste0("O=",donor,"_merged.bam"))
system(picard.code)
