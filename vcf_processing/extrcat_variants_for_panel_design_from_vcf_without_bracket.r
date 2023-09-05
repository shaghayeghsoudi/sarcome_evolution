#rm(list = ls())

####################################################################################
### extrcat variants for the panel design (from vcf files genearted by survivor) ###
####################################################################################

#library(webr)
library(dplyr)
#library(vcfR)
#library(tidyverse)
#library(StructuralVariantAnnotation)
library(VariantAnnotation)
#library(data.table)
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
        mutate(sample_info=gsub("ways","caller",sample_info)) %>% 
        dplyr::select(-raw_id)
       


both<-mutate(type_data_inforow,good_mat_fin) %>%     ### bind two dataframes and remove duplicated columns
      mutate(sample=gsub("_.*$","",sample_info)) %>% 
      mutate(sample=gsub("1","GCT",
                    gsub("2","RD",
                    gsub("3","SW",sample))))

   

chroms<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")   
both<-both[(both$seqnames%in%chroms) & (both$CHR2%in%chroms),]
#both$caller_count_info<-sub('.*\\_', '', both$sample_info)


samples<-unique(both$sample)
callers<-c("1","2","3","4")
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

        #focal_caller<-focal[(focal$caller_count_info==paste(jj,"caller",sep ="") & focal$caller_count==jj),]  ### caller_count column
        focal_caller<-focal[focal$SUPP==callers[jj],]
        
        out_res_type<-NULL
        for (ii in 1:length(var_types)) { ### loop through each variant type

            focal_caller_vartype<-focal_caller[focal_caller$SVTYPE==var_types[ii],]
            if (var_types[ii]=="INV"){
        
            focal_caller_vartype<-focal_caller_vartype[focal_caller_vartype$nanoSV!=1,] ### nanoSV inverions are not inversion
             }

             #if (var_types[ii]=="INS"){ 
             #
            # focal_caller_vartype$END<-(focal_caller_vartype$END)
             
            #}

            #if ((nrow(focal_caller_vartype)) >=57) {
            #
            #focal_random<-focal_caller_vartype[sample(nrow(focal_caller_vartype), 57+15,replace = TRUE) ,]
#
             #} 
        
           #if ((nrow(focal_caller_vartype)) < 57) {

            #focal_random<-focal_caller_vartype[sample(nrow(focal_caller_vartype),nrow(focal_caller_vartype), replace = TRUE) ,]

            # }    


           focal_random<-focal_caller_vartype
            
            out_res_type<-rbind(focal_random,out_res_type) 
            #out_res_type<-out_res_type[!duplicated(out_res_type$start_end),]
            #out_res_type<-out_res_type[,c("seqnames","start", "end","width","strand", "paramRangeID","REF", "ALT","QUAL","FILTER" ,"sample_info", "CIEND", "CIPOS" , "CHR2" ,"END" ,"MAPQ" ,"RE", "IMPRECISE" ,"PRECISE" ,"SVLEN" ,"SVMETHOD","SVTYPE","SUPP_VEC","SUPP","STRANDS","CuteSV", "nanoSV" ,"Sniffles" ,"Svim" ,"sample_info.1", "caller_count" ,"sample" ,"start_end" ,"f" ,"r")]

           } ### variant type loop

            #if(nrow(out_res_type)<284) {

            #additions<-focal_caller[focal_caller$SVTYPE=="DEL",]
            #additions_good<-additions[!(additions$start%in%out_res_type$start),]
            #additions_random<-additions_good[sample(nrow(additions_good), 284-nrow(out_res_type),replace = TRUE) ,]
            #additions_random<-additions_good[sample(nrow(additions_good), 330-nrow(out_res_type),replace = TRUE) ,]

            #out_res_type<-rbind(out_res_type,additions_random)

        #} 
        

    out_res<-rbind(out_res_type,out_res)

    } ## number of callers loop

    #out_res[out_res$SVTYPE=="INS","END"]<-out_res[out_res$SVTYPE=="INS","start"]+out_res[out_res$SVTYPE=="INS","SVLEN"]
    out_res_final<-out_res[,c("SVTYPE","seqnames","start","f","CHR2","END","r","ALT","CuteSV", "nanoSV", "Sniffles", "Svim" )]
    out_res_final_dedup<-out_res_final %>% distinct()

    setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels")
    write.table(out_res_final_dedup,file = paste(samples[kk],"pacbioSV.bed",sep = "_"),col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)



}  ## sample loop
    
###########
### END ###
###########