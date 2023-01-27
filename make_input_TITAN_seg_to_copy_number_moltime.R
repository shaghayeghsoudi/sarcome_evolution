#!/usr/bin/env Rscript

### This script makes copy_number.txt file for mol_time analysis from TITAN Seg file

#rm(list = ls())
#library(tidyverse)
##### ATTRIBUTION #####
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 


### example input from mol_time documenttaion
#sample	Chrom	start	end	class	nMaj1_A	nMin1_A	frac1_A	nMaj2_A	nMin2_A	code_batt
#PD26400a	1	762601	84710840	normal	1	1	1	NA	NA	no_sub
#PD26400a	1	84714419	94909862	deletion	1	0	1	NA	NA	no_sub

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/")
### list and load input files 
filenames <- list.files("SegsFiles_test",pattern="*.segs.txt", full.names = TRUE)
attackStats <- lapply(filenames,function(x) {
      read.csv(x,  header=TRUE, sep = "\t")[,c(1:4,8:14)]
    })

### add sample name as a column
#for (i in 1:length(attackStats)){
#    attackStats[[i]]<-cbind(attackStats[[i]],filenames[i])
#   }
aa <- do.call("rbind", attackStats) 
#head(aa)

aa$Chromosome<-gsub("chr","",aa$Chromosome)
#aa$zz<-gsub("SegsFiles_test/", "",gsub("_cluster2.segs.txt","",aa[,12]))

### make a column to assign clonal and sub_clonal status
aa$code_batt = ifelse(aa$Cellular_Prevalence==1,"no_sub","subclonal")

### make a column to assign normal, gain, deletion, LOH, WCD status
aa$class<- ifelse(aa$code_batt== "no_sub", 
    ifelse (aa$MajorCN==1 & aa$MinorCN==1, "normal",
    ifelse(aa$MajorCN<=1 & aa$MinorCN==0, "deletion",
    ifelse(aa$MajorCN==2 & aa$MinorCN==1, "gain",
    ifelse(aa$MajorCN==3 & aa$MinorCN==1, "gain_2",
    ifelse(aa$MajorCN==4 & aa$MinorCN==1, "gain_3", 
    ifelse(aa$MajorCN==5 & aa$MinorCN==1, "gain_high",
    ifelse(aa$MajorCN>=2 & aa$MinorCN>=2, "Whole_chromosome_duplication",
    ifelse(aa$MajorCN==2 & aa$MinorCN==0, "LOH",
    ifelse(aa$MajorCN==3 & aa$MinorCN==0, "LOH_2",
    ifelse(aa$MajorCN==4 & aa$MinorCN==0, "LOH_3",
    ifelse(aa$MajorCN>=5 & aa$MinorCN==0, "LOH_High","-")
    )))))))))),"-")
 


## adjust colnames and order colmns according to the maunal
aa_fix<-aa[!(aa$MajorCN==0 & aa$MinorCN==0),] 
aa_clonal<-aa_fix[aa_fix$class!="subclonal",]

aa_clonal$nMaj2_A<-NA
aa_clonal$nMin2_A<-NA

needed_columns<-aa_clonal[,c("Sample","Chromosome","Start_Position.bp." ,"End_Position.bp.","class","MinorCN", "MajorCN","Cellular_Prevalence","nMaj2_A","nMin2_A","code_batt")]
colnames(needed_columns)<-c("sample","Chrom","start","end","class","nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","code_batt")
ind <- apply(needed_columns, 1, function(x) all(is.na(x)))
needed_columns_noNA <- needed_columns[ !ind, ]
needed_columns_noNA$sample<-gsub("-.*$", "", needed_columns_noNA$sample)   ## to make sample names consistent with sample names in the ploidy file



write.table(needed_columns_noNA, file = "moltime_data/data/copy_number.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")

#### format ploidy-purity file

ploidy<-read.table(file = "Titan_allsamples_purity_ploidy_test.txt", header = FALSE)[,c(1,2,3)]
colnames(ploidy)<-c("SAMPLE","ABBR_CELL_FRAC","TUM_PLOIDY")
write.table(ploidy, file = "moltime_data/data/sample_normal_contamination.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")

## running command
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")
mol_time(data.dir="moltime_data/data", res.dir="moltime_data/res-defaults")