
library(GenomicRanges)

#gistic2 -b ./softwares/test/output_test/ -mk ./softwares/test/examplefiles/markersfile.txt -seg ./softwares/test/examplefiles/segmentationfile.txt -refgene ./softwares/gistic2/support/refgenefiles/hg19.mat 
#gistic2 -b ./softwares/test/output_sarcoma  -seg ./softwares/test/example_sarcoma/cnv_data_sarcoma_segmenttaion.txt -refgene ./softwares/gistic2/support/refgenefiles/hg19.mat 
gistic2 -b ./99-base_outputs-merged-per-sample  -seg ./03-inputs-gistic/cnv_data_sarcoma_merged_across_regions_meanCNA.txt -refgene /home/users/shsoudi/softwares/gistic2/support/refgenefiles/hg19.mat 
gistic2 -b ./03-inputs-gistic/test2_out  -seg ./03-inputs-gistic/test2/ditest2_ordered_final.txt -refgene /home/users/shsoudi/softwares/gistic2/support/refgenefiles/hg19.mat 


#### prepare input segmentation file for running GISTIC2
setwd("/scratch/users/shsoudi/analysis/gistic")
list_cnvs<-list.files ("00-inputs-segs", pattern = "*.segs.txt", full.names = TRUE)
seg_files<-lapply(list_cnvs, function(x){
    read.delim(x, header=TRUE, sep = "\t")[,c(1,2,3,4,5,10)]
})


cnv_data <- do.call("rbind", seg_files) 
cnv_data$Sample<-gsub("-.*$","",cnv_data$Sample)
cnv_data$Chromosome<-gsub("chr","",cnv_data$Chromosome)
colnames(cnv_data)[c(3:4)]<-c("Start","End")
cnv_data$Seg.CN<- (log2(cnv_data$Copy_Number) -1 )

cnv_data<-cnv_data[,c(1,2,3,4,5,7)]
write.table(cnv_data, file = "01-inputs-gistic/cnv_data_sarcoma_all_segmenttaion.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)



meta<-read.table(file = "/home/users/shsoudi/metadata_updated.txt", header = FALSE)

cnv_data_RTid<-merge(cnv_data,meta, by.x ="Sample" , by.y = "V1")
cnv_in_meta<-cnv_data_RTid[,c(1:6)]
write.table(cnv_in_meta, file = "01-inputs-gistic/cnv_data_sarcoma_all_in_metadata_segmenttaion.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


postRT<-cnv_data_RTid[grepl("postRT",cnv_data_RTid$V5),c(1:6)]
write.table(postRT, file = "01-inputs-gistic/cnv_data_sarcoma_postRT_segmenttaion.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

noRT<-cnv_data_RTid[!grepl("postRT",cnv_data_RTid$V5),c(1:6)]
write.table(noRT, file = "01-inputs-gistic/cnv_data_sarcoma_noRT_segmenttaion.txt", col.names  = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


#####
maf_files<-read.delim(file = "mafs/all_samples_updated_converted_like.maf", sep = "\t")
samples<-unique(cnv_data$Sample)


for (ii in 1:length(samples)){

    genomic_ranges<-cnv_data[cnv_data$Sample==samples[ii],]
    cols <- c("Chromosome","Start", "End" )
    genomic_ranges$id <- apply( genomic_ranges[ , cols ] , 1 , paste , collapse = "-" )
    snp_table<-maf_files[maf_files$Tumor_Sample_Barcode ==samples[ii],]
    cnv_gr<-GRanges(IRanges(start = qq$Start,end =  qq$End, names =qq$id ),seqnames = qq$Chromosome)
    pps_gr<-GRanges(IRanges(start = as.numeric(snp_table$Start_Position),end =  as.numeric(snp_table$End_Position)),seqnames = snp_table$Chromosome)
    count_per_segment<-data.frame("Num Markers"=countOverlaps( cnv_gr,pps_gr))
    count_per_segment$seg_id<-rownames(count_per_segment)
    rownames(count_per_segment)<-NULL
    final<-merge(qq,count_per_segment, by.x = "id", by.y ="seg_id" )
    final$Seg.CN<- (log2(final$Copy_Number) -1 )


}





cnv_gr<-GRanges(IRanges(start = qq$Start,end =  qq$End, names =qq$id ),seqnames = qq$Chromosome)
pps_gr<-GRanges(IRanges(start = as.numeric(snp_table$Start_Position),end =  as.numeric(snp_table$End_Position)),seqnames = snp_table$Chromosome)

countOverlaps(gr, grl)






