#!/usr/bin/env Rscript

#rm(list = ls())
## this script convert maf file into require SNV file format (i.e., shared_clonal_muttaions.txt) required by mol_time package
library(tidyverse)

##### ATTRIBUTION #####
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 


### convert maf to mol_time ###

### example input from mol_time documentation (shared_clonal_muttaion.txt)
#Sample	Chrom	Pos	Ref	Alt	Gene	PM.Tum
#PD26400a	1	102743557	C	T	RP11-202K23.1	1

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")
### load and read all maf files
filenames_CN <- list.files("filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.csv(x, header=TRUE, sep = "\t")[,c("sample_id","Chromosome","Start","Ref", "Alt","Gene","AF")]
     })

## calulate VAF
#zz<-lapply(attackStats_CN, function(x) transform(x, PM.Tum = (t_alt_count/(t_alt_count+t_ref_count))))

### add a column as sample name
#for (i in 1:length(zz)){
#    zz[[i]]<-cbind(zz[[i]],filenames_CN[i])
#    }

shared_clonal_mutations <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix<-shared_clonal_mutations[!rowSums(nchar(as.matrix(shared_clonal_mutations[4:5]))!=1),]
shared_clonal_mutations_fix<-shared_clonal_mutations_fix[!grepl("-", shared_clonal_mutations_fix$Reference_Allele) & !grepl("-", shared_clonal_mutations_fix$Tumor_Seq_Allele2),]


#shared_clonal_mutations_fix$Sample<-gsub("genome-hg19/", "",gsub("--matched.final.maf","",shared_clonal_mutations_fix[,9]))
colnames(shared_clonal_mutations_fix)<-c("Sample","Chrom","Pos","Ref","Alt","Gene","PM.Tum")
shared_clonal_mutations_fix$Chrom<-gsub("chr","",shared_clonal_mutations_fix$Chrom)
shared_clonal_mutations_fix<-shared_clonal_mutations_fix%>%drop_na()
shared_clonal_mutations_fix<-shared_clonal_mutations_fix[!grepl("-", shared_clonal_mutations_fix$Ref) & !grepl("-", shared_clonal_mutations_fix$Alt),]


write.table(shared_clonal_mutations_fix, file = "moltime/shared_clonal_mutations.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")
