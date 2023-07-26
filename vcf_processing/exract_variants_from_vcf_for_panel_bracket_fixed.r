#rm(list = ls())

### extrcat variants for the panel design (from vcf files genearted by survivor)###

###################################################
##### analyze survivor consensus vcf files ########
###################################################
## load required libraries
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


### load vcf survivor files
files_survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor", pattern = "*.vcf", full.names = TRUE,recursive=TRUE)

### load rowranges field
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

#for (i in 1:length(vcfs_info)){
#    vcfs_info[[i]]<-cbind(vcfs_info[[i]],files_survivor[i])
#    }


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
both$caller_count_info<-sub('.*\\_', '', both$sample_info)


samples<-unique(both$sample)
callers<-c("1","2","3")
var_types<-c(unique(both$SVTYPE))


for (kk in 1:length(samples)){  ### loop through each sample

    focal<-both[both$sample==samples[kk],] %>% 
    mutate(start_end=paste(start, start+SVLEN, sep = "-")) 
    #filter(!duplicated(start_end)) 

    focal$f<-substr(focal$STRANDS, 1,1)
    focal$r<-substr(focal$STRANDS, 2,2)
    #focal_good<-focal[focal$sample_info=="1_SVLEN100_DIS100_1caller",]  ### take everything from 1caller as it has all combinations in it

    
    out_res<-NULL
    for(jj in 1:length(callers)){  ### loop through each number of variant caller detected a variant (caller_count_info column)

        focal_caller<-focal[(focal$caller_count_info==paste(jj,"caller",sep ="") & focal$caller_count==jj),]  ### caller_count column

        
        out_res_type<-NULL
        for (ii in 1:length(var_types)) { ### loop through each variant type

            focal_caller_vartype<-focal_caller[focal_caller$SVTYPE==var_types[ii],]
            if (var_types[ii]=="INV"){
        
            focal_caller_vartype<-focal_caller_vartype[focal_caller_vartype$nanoSV!=1,]
             }

            if ((nrow(focal_caller_vartype)) >=57) {

            focal_random<-focal_caller_vartype[sample(nrow(focal_caller_vartype), 57+15,replace = TRUE) ,]

             } 
        
           if ((nrow(focal_caller_vartype)) < 57) {

            focal_random<-focal_caller_vartype[sample(nrow(focal_caller_vartype),nrow(focal_caller_vartype), replace = TRUE) ,]

             }    

            out_res_type<-rbind(focal_random,out_res_type) 
            #out_res_type<-out_res_type[!duplicated(out_res_type$start_end),]
        #out_res_type<-out_res_type[,c("seqnames","start", "end","width","strand", "paramRangeID","REF", "ALT","QUAL","FILTER" ,"sample_info", "CIEND", "CIPOS" , "CHR2" ,"END" ,"MAPQ" ,"RE", "IMPRECISE" ,"PRECISE" ,"SVLEN" ,"SVMETHOD","SVTYPE","SUPP_VEC","SUPP","STRANDS","CuteSV", "nanoSV" ,"Sniffles" ,"Svim" ,"sample_info.1", "caller_count" ,"sample" ,"start_end" ,"f" ,"r")]

        } ### variant type loop

        if(nrow(out_res_type)<284) {

            additions<-focal_caller[focal_caller$SVTYPE=="DEL",]
            additions_good<-additions[!(additions$start%in%out_res_type$start),]
            #additions_random<-additions_good[sample(nrow(additions_good), 284-nrow(out_res_type),replace = TRUE) ,]
            additions_random<-additions_good[sample(nrow(additions_good), 330-nrow(out_res_type),replace = TRUE) ,]

            out_res_type<-rbind(out_res_type,additions_random)

        } 
        

    out_res<-rbind(out_res_type,out_res)

    } ## number of callers loop

    #out_res[out_res$SVTYPE=="INS","END"]<-out_res[out_res$SVTYPE=="INS","start"]+out_res[out_res$SVTYPE=="INS","SVLEN"]
    out_res_final<-out_res[,c("SVTYPE","seqnames","start","f","CHR2","END","r","ALT","CuteSV", "nanoSV", "Sniffles", "Svim" )] ## ALT column included for bracket information
    out_res_final_dedup<-out_res_final %>% distinct()

    setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels")
    write.table(out_res_final_dedup,file = paste(samples[kk],"pacbioSV.bed",sep = "_"),col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)



}  ## sample loop
    
###################
###################
## To fix brackets

chroms<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")   


### load original vcf files
files_original<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels/filtered_DV_5_vcfs/site_annotations", pattern = "*.txt", full.names = TRUE,recursive=TRUE)

### load rowranges field
attackStats <- lapply(files_original,function(x) {
      read.csv(x,  header=FALSE, sep = "\t")
    })

for (i in 1:length(attackStats )){
    attackStats [[i]]<-cbind(attackStats [[i]],files_original[i])
    }
type_data <- do.call("rbind", attackStats) 
type_data$SVTYPE<-"TRA"

names(type_data)[4]<-"path"
type_data_good<-type_data %>% mutate(sample_info=sub('.*/\\s*', '',gsub(".vcf.info.txt","",path))) %>% 
      mutate(caller_info = gsub(".*_","",sample_info)) %>% 
      mutate(sample_info_good = gsub("_.*$","",sample_info)) %>% 
      filter(V1%in%chroms)

type_data_good$sample_info_good<-ifelse(type_data_good$sample_info_good==1, "GCT",
                                 ifelse(type_data_good$sample_info_good==2, "RD",
                                 ifelse(type_data_good$sample_info_good==3, "SW","-")))
                                 
                                  


### read bed like files 
files_original_bed<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels", pattern = "*.txt", full.names = TRUE)

attackStats_bed <- lapply(files_original_bed,function(x) {
      read.csv(x,  header=TRUE, sep = "\t")
    })

for (i in 1:length(attackStats_bed )){
    attackStats_bed [[i]]<-cbind(attackStats_bed [[i]],files_original_bed[i])
    }
bed_data <- do.call("rbind", attackStats_bed) 

names(bed_data )[13]<-"path"
bed_data_good<-bed_data %>% 
     mutate(sample_info=sub('.*/\\s*', '',gsub("_pacbioSV.with.bracetInfo.072523.txt","",path))) %>% 
     filter(seqnames%in%chroms)


#bed_data_tra<-bed_data_good[bed_data$SVTYPE=="TRA",]

samples<-unique(bed_data_tra$sample_info)

### add fake ID for the bed file
out_res_bed<-NULL
for (ii in 1:length(samples)){

    ### process bed file
    bed_data_good_focal<-bed_data_good[bed_data_good$sample_info==samples[ii],] 
    bed_data_good_focal_tra<-bed_data_good_focal[bed_data_good_focal$SVTYPE=="TRA",]
    bed_data_good_focal_tra_fake<-bed_data_good_focal_tra %>% 
    mutate(start_fake=substr(start,1,nchar(start)-1)) %>%   ### droping the last digit
    mutate(end_fake=substr(END,1,nchar(END)-1)) %>% 

    filter(CHR2%in%chroms) %>%
    mutate(fake_id= paste(seqnames,start_fake,CHR2,end_fake, sep = "_")) %>% 
    rename(CuteSV = "caller_CuteSV", nanoSV= "caller_nanoSV", Sniffles="caller_Sniffles" , Svim = "caller_Svim") 

    tra_fake<-bed_data_good_focal_tra_fake%>% mutate(min_wtp = apply(bed_data_good_focal_tra_fake[, grepl("caller", names(bed_data_good_focal_tra_fake))], 1, function(x) {
    names(x)[min(which(x > 0))]  ### Rowwise name of column where first non-zero value appears
      }))

      out_res_bed<-rbind(tra_fake,out_res_bed)
}
out_res_bed_good<-dplyr::select(out_res_bed,-path) %>% 
    
    
    

out_res_type<-NULL
for (ii in 1:length(samples)){   
    ### type data
    type_data_good_focal<-type_data_good[type_data_good$sample_info_good==samples[ii],]  ### find focal sample name
    #type_data_good_focal<-type_data_good_focal[type_data_good_focal$caller_info=="cuteSV",]
    type_data_good_focal_fake<-type_data_good_focal %>% 
    mutate(CHR2=gsub(":.*$","",V3),end =gsub(".*:","",V3)) %>% 
    mutate(CHR2 = str_split(CHR2, 'chr', simplify = TRUE)[,2]) %>% 
    mutate(end = as.numeric(gsub("\\D+", "", end))) %>% 
    mutate(CHR2= paste("chr",CHR2,sep = "")) %>% 
    filter(CHR2%in%chroms, V1%in%chroms) %>%
    mutate(end_fake=substr(end,1,nchar(end)-1),V2=as.numeric(V2)) %>% 
    mutate(start_fake=as.numeric(substr(V2,1,nchar(as.numeric(V2))-1))) %>% 
    mutate(fake_id= paste(V1,start_fake,CHR2,end_fake, sep = "_"))
    out_res_type<-rbind(type_data_good_focal_fake,out_res_type)

}
out_res_type_good<-dplyr::select(out_res_type,-path) 
     


### assigns the right bracket from original vcf cite annotation
for (kk in 1:length(samples)){


    out_res_bed_good_focal<-out_res_bed_good[out_res_bed_good$sample_info==samples[kk],]  ### focal sample bed like format
    out_res_bed_good_focal$min_wtp<- gsub("caller_CuteSV", "cuteSV",
           gsub("caller_nanoSV", "nanoSV",
           gsub("caller_Sniffles", "sniffles", 
           gsub("caller_Svim","svim",out_res_bed_good_focal$min_wtp))))


    out_res_type_good_focal<-out_res_type_good[out_res_type_good$sample_info_good==samples[kk],]
   

    out_res_pp<-NULL
    for (pp in 1:nrow(out_res_bed_good_focal)){

        out_res_bed_good_focal_line<-out_res_bed_good_focal[pp,] 
        out_res_type_good_focal[out_res_type_good_focal$caller_info==out_res_bed_good_focal_line$min_wtp,]

        ### find overlapping rows with bed file
        qq<-out_res_type_good_focal[out_res_type_good_focal$caller_info==out_res_bed_good_focal_line$min_wtp & out_res_type_good_focal$fake_id==out_res_bed_good_focal_line$fake_id,]
        out_res_pp<-rbind(qq,out_res_pp)

    }


}





##############
#### END #####
##############

