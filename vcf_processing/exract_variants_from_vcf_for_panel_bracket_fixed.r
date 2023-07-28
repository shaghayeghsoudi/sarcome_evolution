#rm(list = ls())

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
files_survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor", 
pattern = "*.vcf", full.names = TRUE,recursive=TRUE)

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
    mutate(sample_info=sub('.*/\\s*', '', gsub("survivor_merged_filtered_PacBio_","",path)), ) 
%>% 
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
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor", 
pattern = "*.txt", full.names = TRUE,recursive=TRUE)


matrix_tables<-lapply(matrix,function(x){
    read.csv(x, header = FALSE, strip.white=T, sep='')
})

for (i in 1:length(matrix_tables)){
    matrix_tables[[i]]<-cbind(matrix_tables[[i]],matrix[i])
    }
mat <- do.call("rbind", matrix_tables) 
names(mat)[5]<-"raw_id"

good_mat<-mat %>% mutate(sample_info=sub('.*/\\s*', 
'',gsub("survivor_merged_filtered_PacBio_","",raw_id)))
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
    #focal_good<-focal[focal$sample_info=="1_SVLEN100_DIS100_1caller",]  ### take everything from 
1caller as it has all combinations in it

    
    out_res<-NULL
    for(jj in 1:length(callers)){  ### loop through each number of variant caller detected a 
variant (caller_count_info column)

        focal_caller<-focal[(focal$caller_count_info==paste(jj,"caller",sep ="") & 
focal$caller_count==jj),]  ### caller_count column

        
        out_res_type<-NULL
        for (ii in 1:length(var_types)) { ### loop through each variant type

            focal_caller_vartype<-focal_caller[focal_caller$SVTYPE==var_types[ii],]
            if (var_types[ii]=="INV"){
        
            focal_caller_vartype<-focal_caller_vartype[focal_caller_vartype$nanoSV!=1,]
             }

             #if (var_types[ii]=="INS"){ 
             #
            # focal_caller_vartype$END<-(focal_caller_vartype$END)
             
            #}

            #if ((nrow(focal_caller_vartype)) >=57) {
            #
            #focal_random<-focal_caller_vartype[sample(nrow(focal_caller_vartype), 57+15,replace 
= TRUE) ,]
#
             #} 
        
           #if ((nrow(focal_caller_vartype)) < 57) {

            
#focal_random<-focal_caller_vartype[sample(nrow(focal_caller_vartype),nrow(focal_caller_vartype), 
replace = TRUE) ,]

            # }    


           focal_random<-focal_caller_vartype
            
            out_res_type<-rbind(focal_random,out_res_type) 
            #out_res_type<-out_res_type[!duplicated(out_res_type$start_end),]
        #out_res_type<-out_res_type[,c("seqnames","start", "end","width","strand", 
"paramRangeID","REF", "ALT","QUAL","FILTER" ,"sample_info", "CIEND", "CIPOS" , "CHR2" ,"END" 
,"MAPQ" ,"RE", "IMPRECISE" ,"PRECISE" ,"SVLEN" 
,"SVMETHOD","SVTYPE","SUPP_VEC","SUPP","STRANDS","CuteSV", "nanoSV" ,"Sniffles" ,"Svim" 
,"sample_info.1", "caller_count" ,"sample" ,"start_end" ,"f" ,"r")]

        } ### variant type loop

        #if(nrow(out_res_type)<284) {

            #additions<-focal_caller[focal_caller$SVTYPE=="DEL",]
            #additions_good<-additions[!(additions$start%in%out_res_type$start),]
            #additions_random<-additions_good[sample(nrow(additions_good), 
284-nrow(out_res_type),replace = TRUE) ,]
            #additions_random<-additions_good[sample(nrow(additions_good), 
330-nrow(out_res_type),replace = TRUE) ,]

            #out_res_type<-rbind(out_res_type,additions_random)

        #} 
        

    out_res<-rbind(out_res_type,out_res)

    } ## number of callers loop

    
#out_res[out_res$SVTYPE=="INS","END"]<-out_res[out_res$SVTYPE=="INS","start"]+out_res[out_res$SVTYPE=="INS","SVLEN"]
    out_res_final<-out_res[,c("SVTYPE","seqnames","start","f","CHR2","END","r","ALT","CuteSV", 
"nanoSV", "Sniffles", "Svim" )]
    out_res_final_dedup<-out_res_final %>% distinct()

    
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels")
    write.table(out_res_final_dedup,file = paste(samples[kk],"pacbioSV.bed",sep = "_"),col.names 
= TRUE, row.names = FALSE, sep = "\t",quote = FALSE)



}  ## sample loop
    

#### END #####


## from out_res
out_res_final_dedup

bracet<-data.frame("brachet_info"=alt<-unlist(out_res_final_dedup$ALT))

out_res_bracet<-cbind(out_res_final_dedup,bracet)
out_res_bracet$identity<-paste(out_res_bracet$SVTYPE ,out_res_bracet$seqnames, 
out_res_bracet$start ,out_res_bracet$f ,out_res_bracet$CHR2 ,out_res_bracet$END, out_res_bracet$r 
,out_res_bracet$CuteSV ,out_res_bracet$nanoSV , out_res_bracet$Sniffles, out_res_bracet$Svim,sep 
= "_")

out_res_bracet_inv<-out_res_bracet[out_res_bracet$SVTYPE=="INV",]
out_res_bracet_TRA<-out_res_bracet[out_res_bracet$SVTYPE=="TRA",]
out_res_bracet_ins<-out_res_bracet[out_res_bracet$SVTYPE=="INS",]


###
sample<-read.delim(file = 
"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels/GCT_pacbioSV.bed", 
header = TRUE)
sample[sample$SVTYPE=="INS","END"]<-sample[sample$SVTYPE=="INS","start"]
sample$identity<-paste(sample$SVTYPE,sample$seqnames,sample$start,sample$f , sample$CHR2, 
sample$END, sample$r,sample$CuteSV,sample$nanoSV,sample$Sniffles,sample$Svim, sep ="_")

sample_inv<-sample[sample$SVTYPE=="INV",]
sample_tra<-sample[sample$SVTYPE=="TRA",]
sample_ins<-sample[sample$SVTYPE=="INS",]


sample_not_invtra<-sample[!(sample$SVTYPE%in%sample_inv$SVTYPE) 
&!(sample$SVTYPE%in%sample_tra$SVTYPE) ,]
sample_not_invtra$brachet_info<-"NA"  ### final all but not tra/inv


###
merge_inv<-merge(sample_inv , out_res_bracet_inv, by.x = "identity", by.y = "identity")
merge_tra<-merge(sample_tra , out_res_bracet_TRA, by.x = "identity", by.y = "identity")
merge_ins<-merge(sample_ins , out_res_bracet_ins, by.x = "identity", by.y = "identity")


merged_inv_tra<-rbind(merge_inv,merge_tra)

merged_inv_tra<-merged_inv_tra[,c("SVTYPE.x","seqnames.x", "start.x","f.x" ,"CHR2.x"   ,  "END.x" 
,"r.x" ,"CuteSV.x", "nanoSV.x" ,"Sniffles.x", "Svim.x" ,"identity" ,"brachet_info")]
names(merged_inv_tra)<-c("SVTYPE","seqnames","start","f","CHR2","END","r","CuteSV", "nanoSV", 
"Sniffles", "Svim","identity", "brachet_info")


final_both<-rbind(sample_not_invtra,merged_inv_tra)
final_both_format<-final_both[,c("SVTYPE","seqnames","start","f","CHR2","END","r","CuteSV", 
"nanoSV", "Sniffles", "Svim","brachet_info")]


write.table(final_both_format, file = 
"Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels/SW_pacbioSV.with.bracetInfo.072523.txt",col.names 
= TRUE, row.names = FALSE, sep = "\t",quote = FALSE)



qq<-read.delim(file = "SW_pacbioSV.with.bracetInfo.072523.txt", header = TRUE)
qq[qq$SVTYPE=="INS","END"]<-qq[qq$SVTYPE=="INS","END"]+1
write.table(qq, file = 
"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor/panels/SW_pacbioSV.with.bracetInfo.072523.txt",col.names 
= TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

