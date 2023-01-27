

############################################################
##### process mutations and run dNdScv fpr Sarcoma samples (all samples) #####
#rm(list = ls())
library("tidyverse")
library("dndscv")
library("plyr")
library("VennDiagram")
library("maftools")

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")

## selected samples
#meta<-read.table(file = "selected_samples.txt", header = FALSE, sep= "\t")
#colnames(meta)<-c("sample_id","purity","ploidy","RT_status","RT_status2","RT_staus_code")
#colnames(meta)<-c("sample_id","RT_status","RT_staus_code")

## all samples
meta<-read.table(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)

## subset to samples that 
meta_selected_prepos<-meta[meta$unique_sample_id%in%c("SRC125","SRC127","SRC167","SRC169","SRC170","SRC171","TB9051"),]

### load and read all maf files
filenames_CN <- list.files("filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.csv(x, header=TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange")]
     })


shared_clonal_mutations <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix<-shared_clonal_mutations%>%drop_na(sample_id)
shared_clonal_mutations_fix$Chromosome<-gsub("chr","",shared_clonal_mutations_fix$Chromosome)

#mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange","RT_status","RT_staus_code")]
mutaion_RT<-merge(shared_clonal_mutations_fix,meta_selected_prepos, by.x = "sample_id", by.y = "sample_id")
mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]


mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

mutaion_RT$RTtretment<-ifelse(mutaion_RT$RT_status== "preRT", "naive",
                        ifelse(mutaion_RT$RT_status== "noRT", "naive",
                        ifelse(mutaion_RT$RT_status== "postRT", "treatment",
                        "-")))
