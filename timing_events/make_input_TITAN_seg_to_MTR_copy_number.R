#!/usr/bin/env Rscript

### This script makes copy_number Granges file for mutationtimeR package from TITAN Seg file

#rm(list = ls())
#library(tidyverse)
library(GenomicRanges)
##### ATTRIBUTION #####
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 

setwd("/data/users/soudi/sarcoma/sarcoma_inputs")
ploidy<-read.delim(file = "Titan_allsamples_purity_ploidy_add500_purity.1_sort.txt",header =FALSE)
colnames(ploidy)<-c("sample","purity","ploidy","sex")


filenames <- list.files("SegsFiles_test",pattern="*.segs.txt", full.names = TRUE)
attackStats <- lapply(filenames,function(x) {
      read.csv(x,  header=TRUE, sep = "\t")[,c("Sample","Chromosome","Start_Position.bp.","End_Position.bp.","MajorCN","MinorCN","Cellular_Prevalence")]
    })


aa <- do.call("rbind", attackStats)    
aa_clonal<-aa[aa$Cellular_Prevalence==1,]
aa_clonal$Chromosome<-gsub("chr","", aa_clonal$Chromosome)
aa_clonal$Sample<-gsub("-.*$", "", aa_clonal$Sample)   ## to make sample names consistent with sample names in the ploidy file

aa_required_cols<-merge(aa_clonal,ploidy, by.x = "Sample", by.y = "sample")[,c(1:6,8)]

ind <- apply(aa_required_cols, 1, function(x) all(is.na(x)))
aa_required_cols_noNA <- aa_required_cols[ !ind, ]


gr_cna_seg <- GRanges(ranges=IRanges(start=aa_required_cols_noNA$Start_Position.bp., end=aa_required_cols_noNA$End_Position.bp.), 
              seqnames=aa_required_cols_noNA$Chromosome,
              strand="*",
              major_cn= aa_required_cols_noNA$MajorCN,
              minor_cn= aa_required_cols_noNA$MinorCN,
              clonal_frequency = aa_required_cols_noNA$purity)