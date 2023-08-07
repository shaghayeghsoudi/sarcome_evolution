#rm(list = ls())


=##################################
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
library("plyranges")
#library(stringi)

#### load hg19 cosmic cancer genes
cosmic<-read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/cancer_genes/Census_May2023_hg19_fromCosmic.csv", header = TRUE)[,c(1,4)]
cosmic_good<-cosmic %>% separate (Genome.Location,c("chrom","start","end")) %>% 
            mutate(chrom = paste("chr",chrom,sep = "")) %>% 
            drop_na() 

cosmic_gr<-makeGRangesFromDataFrame(cosmic_good, keep.extra.columns=TRUE)


### laod survivor merged vcfs ###
### => TO DO: drop contigs 
files_survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/circus_plot/cell_lines/pacbio_survivor_3matched_SVLEN100_DIS100", pattern = "*.vcf", full.names = TRUE,recursive=TRUE)


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
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/circus_plot/cell_lines/pacbio_survivor_3matched_SVLEN100_DIS100", pattern = "*.txt", full.names = TRUE,recursive=TRUE)

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

both_three_ways<-both %>% 
      mutate(caller_count_info<-sub('.*\\_', '', both$sample_info)) %>% 
      filter(caller_count >= 3) %>%     #### take variants detected by at least three callers
      dplyr::select(seqnames,start, end ,CHR2, END ,SVLEN ,SVTYPE, SUPP ,STRANDS ,CuteSV,nanoSV ,Sniffles ,Svim ,sample_info.1 ,caller_count ,sample) %>% 
      mutate(uniq_id = paste(seqnames, start,end,CHR2, END ,SVLEN ,SVTYPE,sample, sep = "_"))


both_three_ways_TRA<-both_three_ways[both_three_ways$SVTYPE=="TRA",]

both_three_ways_sel<-both_three_ways[,c("seqnames","start","end","sample" , "uniq_id")]
both_three_ways_gr<-makeGRangesFromDataFrame(both_three_ways_sel, keep.extra.columns=TRUE)
sv_cosmic_overlap<-data.frame(join_overlap_inner(both_three_ways_gr,cosmic_gr))   ### subset to variants overlapping with cosmic genes


write.table(sv_cosmic_overlap, file = "sv_cosmic_overlap.txt", col.names = TRUE, row.names = FALSE)

###### find overlap of main SV set with those overlaping with cosmic genes
both_three_ways_cosmic<-both_three_ways[both_three_ways$uniq_id%in%sv_cosmic_overlap$uniq_id,]


both_three_ways_cosmic_tra<-rbind(both_three_ways_cosmic,both_three_ways_TRA)
###########################
##### start pltting ######

samples<-unique(both_three_ways_cosmic_tra$sample)

for (kk in 1:length(samples)){  ### loop through each sample

    focal<-both_three_ways_cosmic_tra[both_three_ways_cosmic_tra$sample==samples[kk],] 

    focal_myc<-focal[,c("seqnames" ,"start" , "end", "SVLEN","SVTYPE")] %>% 
        rename(seqnames = "chrom_myc", start= "start_myc", end = "end_myc")

    focal_myc$color<-ifelse(focal_myc$SVTYPE== "INS", 2 ,
                     ifelse(focal_myc$SVTYPE== "DEL", 3 ,
                     ifelse(focal_myc$SVTYPE== "DUP", 4 ,
                     ifelse(focal_myc$SVTYPE== "INV", 5 ,
                     ifelse(focal_myc$SVTYPE== "TRA", 6 ,
                     "_")))))

                     
    
    focal_nonmyc<-focal[,c("CHR2" ,"END")] %>% 
        mutate(end= rep(END)) %>% 
        rename(CHR2 = "chrom_nonmyc", END= "start_nonmyc", end = "end_nonmyc")


    ### make ordered chromosmes
    chroms_myc<-unique(focal_myc$chrom_myc) 
    chroms_nonmyc<-as.factor(unique(focal_nonmyc$chrom_nonmyc))
    chroms<-data.frame("chr"=unlist(list(chroms_myc,chroms_nonmyc))) %>% 
    mutate(chr = gsub("chr","",chr)) %>% 
    distinct (chr) %>% 
    #arrange(as.numeric(chr)) %>%
    arrange((chr)) %>%
    mutate(chr= paste("chr",chr, sep= ""))


pdf(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/circus_plot/cell_lines/circlize_plot_overlaped_with_Cosmic",samples[kk],".pdf", sep = ""),width = 5, height = 5)
    circos.par(gap.degree = 2, start.degree = 170)
    circos.initializeWithIdeogram(
    plotType = c("ideogram", "labels"), species = "hg19",
    chromosome.index = c(chroms$chr))
    #sector.width = c(3, 1, 1, 1,1))
 

    circos.genomicLink(focal_myc,
                     focal_nonmyc,
                     lwd=3,
                     col = alpha(focal_myc$color,0.8))
                     #col = rand_color(unique(focal_myc$SVTYPE)), transparency = 0.3)

dev.off()

}



################
##### END ######
################



