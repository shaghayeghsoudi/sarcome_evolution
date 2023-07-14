
### make GISTIC input from seg files ###

#library(GenomicRanges)
library(tidyr)
library(tidyverse)
library(dplyr)
library(plyr)
library("maftools")


### assign relevant directory and create a new directory for the final file 
mainDir <- "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/Correct_CNA_segs_newsolutionsfromanish"
subDir <- "01-inputs-gistic"
if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
} else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    
}

###################################
###### load updated metadata ######
###################################
meta_raw<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")

meta_raw$RTstatus<-ifelse(meta_raw$sequenceofsamplevRT== "beforeRT", "noRT",
                        ifelse(meta_raw$sequenceofsamplevRT== "afterRT", "postRT",
                        ifelse(meta_raw$sequenceofsamplevRT== "nopreopRT", "noRT",
                        "-")))

meta<-meta_raw%>%
  #rename(sampleid="meta_id") %>%
  mutate(unique_sample_id=gsub("_.*$","",sampleid)) %>%
  mutate(identifier=paste(unique_sample_id,RTstatus, sep = "_"))



#### prepare input segmentation file for running GISTIC2 ####
### Load TITAN seg files
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/Correct_CNA_segs_newsolutionsfromanish")
list_cnvs<-list.files ("00-inputs-segs", pattern = "*.seg.txt", full.names = TRUE)
seg_files<-lapply(list_cnvs, function(x){
    read.delim(x, header=TRUE, sep = "\t")[,c(1,2,3,4,7,20)]  ## (Sample,Chromosome,Start,End,Length.snp.,Corrected_Copy_Number)
})

cnv_data <- do.call("rbind", seg_files) %>% 
   mutate(Sample=gsub("-.*$","",Sample), Chromosome= gsub("chr","",Chromosome),Seg.CN = (log2(Corrected_Copy_Number) -1 )) %>% 
   select(Sample,Chromosome ,Start, End, Length.snp.,Seg.CN) %>% 
   filter(Sample%in%meta_raw$sampleid) %>% 
   mutate(Chromosome=gsub("X","23",Chromosome))



write.table(cnv_data, file = "cnv_data_sarcoma_all_regions_fromupdated.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

#############
#### END ###
############


#cnv_data_RTid<-merge(cnv_data,meta, by.x ="Sample" , by.y = "V1")
#cnv_in_meta<-cnv_data_RTid[,c(1:6)]
#write.table(cnv_in_meta, file = "01-inputs-gistic/cnv_data_sarcoma_all_in_metadata_segmenttaion.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#
#
#postRT<-cnv_data_RTid[grepl("postRT",cnv_data_RTid$V5),c(1:6)]
#write.table(postRT, file = "01-inputs-gistic/cnv_data_sarcoma_postRT_segmenttaion.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#
#noRT<-cnv_data_RTid[!grepl("postRT",cnv_data_RTid$V5),c(1:6)]
#write.table(noRT, file = "01-inputs-gistic/cnv_data_sarcoma_noRT_segmenttaion.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


#####
#maf_files<-read.delim(file = "mafs/all_samples_updated_converted_like.maf", sep = "\t")
#samples<-unique(cnv_data$Sample)


#for (ii in 1:length(samples)){
#
#    genomic_ranges<-cnv_data[cnv_data$Sample==samples[ii],]
#    cols <- c("Chromosome","Start", "End" )
#    genomic_ranges$id <- apply( genomic_ranges[ , cols ] , 1 , paste , collapse = "-" )
#    snp_table<-maf_files[maf_files$Tumor_Sample_Barcode ==samples[ii],]
#    cnv_gr<-GRanges(IRanges(start = qq$Start,end =  qq$End, names =qq$id ),seqnames = qq$Chromosome)
#    pps_gr<-GRanges(IRanges(start = as.numeric(snp_table$Start_Position),end =  as.numeric(snp_table$End_Position)),seqnames = snp_table$Chromosome)
#   count_per_segment<-data.frame("Num Markers"=countOverlaps( cnv_gr,pps_gr))
#    count_per_segment$seg_id<-rownames(count_per_segment)
#    rownames(count_per_segment)<-NULL
#    final<-merge(qq,count_per_segment, by.x = "id", by.y ="seg_id" )
#    final$Seg.CN<- (log2(final$Copy_Number) -1 )
#
#
#}
#
#
#cnv_gr<-GRanges(IRanges(start = qq$Start,end =  qq$End, names =qq$id ),seqnames = qq$Chromosome)
#pps_gr<-GRanges(IRanges(start = as.numeric(snp_table$Start_Position),end =  as.numeric(snp_table$End_Position)),seqnames = snp_table$Chromosome)
#
#countOverlaps(gr, grl)

### GISTIC running command ###
gistic2 -b ./99-outputs  -seg ./01-input-GISTIC/cnv_data_sarcoma_all_regions_fromupdated.txt -refgene /home/users/shsoudi/softwares/gistic2/support/refgenefiles/hg19.mat -rx {0}
#gistic2 -b ./03-inputs-gistic/test2_out  -seg ./03-inputs-gistic/test2/ditest2_ordered_final.txt -refgene /home/users/shsoudi/softwares/gistic2/support/refgenefiles/hg19.mat 

