### extrcat variants for the panel design ###

#https://krassowski.github.io/complex-upset/articles/Examples_R.html#fill-the-bars
###################################
##### analyze survivor consensus vcf files
###################################
library(ggplot2)
library(webr)
library(dplyr)
library(plyr)
library(vcfR)
library(ggplot2)
library(tidyverse)
library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(UpSetR)
library(ComplexUpset)
library(data.table)
#library(stringi)


### => TO DO: drop contigs 
### load vcf files
files_survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor", pattern = "*.vcf", full.names = TRUE,recursive=TRUE)

### lead rowranges field
vcfs_survivor<-lapply(files_survivor,function(x){
     vcf<-readVcf(x)
     #info(readVcf(x))
     vcf_row<-data.frame(rowRanges(vcf))
     #data.frame(info(readVcf(x)))
})

for (i in 1:length(vcfs_survivor)){
    vcfs_survivor[[i]]<-cbind(vcfs_survivor[[i]],files_survivor[i])
    }
type_data <- do.call("rbind", vcfs_survivor) 
names(type_data)[11]<-"path"

type_data_good<-type_data%>% 
    mutate(sample_info=sub('.*/\\s*', '', gsub("survivor_merged_filtered_PacBio_","",path)), ) %>% 
    mutate(sample_info=gsub(".vcf","",sample_info)) 


#### load vcf info #####
vcfs_info<-lapply(files_survivor,function(x){
     vcf<-readVcf(x)
     #info(readVcf(x))
     #vcf_row<-data.frame(rowRanges(vcf))
     data.frame(info(readVcf(x)))
})

type_data_info <- do.call("rbind", vcfs_info) 
rownames(type_data_info )<-NULL

type_data_inforow<-cbind(type_data_good,type_data_info)  ### final merged set
type_data_inforow<-type_data_inforow[,colnames(type_data_inforow)!="path"] %>% 
    mutate(sample_info=gsub("2callers","2caller",sample_info))

#########################
### load matrix files ###
#########################
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor", pattern = "*.txt", full.names = TRUE,recursive=TRUE)


matrix_tables<-lapply(matrix,function(x){
    read.csv(x, header = FALSE, strip.white=T, sep='')
})

for (i in 1:length(matrix_tables)){
    matrix_tables[[i]]<-cbind(matrix_tables[[i]],matrix[i])
    }
mat <- do.call("rbind", matrix_tables) 
names(mat)[5]<-"raw_id"

good_mat<-mat %>% mutate(sample_info=sub('.*/\\s*', '',gsub("survivor_merged_filtered_PacBio_","",raw_id)))
good_mat$sample_info<-gsub("overlapped_","",gsub("_matrix.txt","",good_mat$sample_info))

good_mat_fin<-good_mat%>%
        rename(V1="CuteSV" , V2="nanoSV" , V3= "Sniffles",V4="Svim") %>% 
        mutate(caller_count=rowSums(.[1:4])) %>% 
        mutate(sample_info=gsub("ways","caller",sample_info)) 
       
good_mat_fin<-good_mat_fin[,colnames(good_mat_fin)!="raw_id"]
       
both<-cbind(type_data_inforow,good_mat_fin) %>% 
      mutate(sample=gsub("_.*$","",sample_info)) %>% 
      mutate(sample=gsub("1","GCT",sample)) %>% 
      mutate(sample=gsub("2","RD",sample)) %>%
      mutate(sample=gsub("3","SW",sample))
   

chroms<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")   
both<-both[both$seqnames%in%chroms,]
both<-both[both$CHR2%in%chroms,]


samples<-unique(both$sample)
callers<-c("1","2","3")

for (kk in 1:length(samples)){

    focal<-both[both$sample==samples[kk],] %>% 
    mutate(start_end=paste(start, start+SVLEN, sep = "-")) %>% 
    filter(!duplicated(start_end)) 

    focal$f<-substr(focal$STRANDS, 1,1)
    focal$r<-substr(focal$STRANDS, 2,2)

    
    out_res<-NULL
    for(jj in 1:length(callers)){

        focal_caller<-focal[focal$caller_count==callers[jj],] 
        pattern1<-c("INV", "DUP", "TRA")
        focal_invdup<-focal_caller[grepl(paste(pattern1, collapse = "|"), focal_caller$SVTYPE),] 

        if ((nrow(focal_invdup)) >=500) {

            focal_invdup_randdom<-focal_invdup[sample(nrow(focal_invdup), 500, replace = FALSE) ,]

        } 
        
        if ((nrow(focal_invdup)) < 500) {

            focal_invdup_randdom<-focal_invdup[sample(nrow(focal_invdup), nrow(focal_invdup), replace = FALSE) ,]
        }    


    
        ## only ins and dels
        focal_delins<-focal_caller[grepl("INS|DEL",focal_caller$SVTYPE),]
        count_dels<-850-nrow(focal_invdup_randdom)
        
        focal_delins_randdom<-focal_delins[sample(nrow(focal_delins), count_dels, replace = FALSE) ,]

        all_types<-rbind(focal_invdup_randdom,focal_delins_randdom)
        out_res<-rbind(all_types,out_res)

    }
    
    out_res_final<-out_res[,c("SVTYPE","seqnames","start","f","CHR2","END","r","CuteSV", "nanoSV", "Sniffles", "Svim" )]
    
    setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels")
    write.table(out_res_final,file = paste(samples[kk],"pacbioSV.bed",sep = "_"),col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

    

}