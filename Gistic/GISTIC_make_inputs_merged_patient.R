
# Title: 
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 
# Date: January 2023

### load required libraries
#library(bedtoolsr)
#library(GenomicRanges)
library(data.table)
library(VariantAnnotation)

## load metadata and TITAN output segment files 
meta<-read.table(file = "/scratch/users/shsoudi/metadata/metadata.txt", header = FALSE)
seg_files<-list.files("00-inputs-segs", pattern = "*.segs.txt", full.names = TRUE)

attackTitan1 <- function(x) {
      data<-read.delim(x,  header=TRUE, sep = "\t")
      data$sample_id<-gsub("-.*$","",data$Sample)
      data<-merge(data,meta, by.x = "sample_id", by.y = "V1")
      data$uniq_id<-gsub("_.*$","",data$sample_id)
      data$uniq_id2<-paste(data$uniq_id,data$V5, sep = "_")
      #data$major_minor_copy<-paste(data$MajorCN,data$MinorCN, sep = "-")
      return(data)
    }

titan_segs <- lapply(seg_files, attackTitan1)
titan_segs <- do.call("rbind", titan_segs)
patients_rt<-unique(titan_segs$uniq_id2)

## load vcfs and write simple output files consisting of only chrom and position for each region 
#vcf_all<-list.files("test", pattern = "*.vcf", recursive = TRUE, full.names = TRUE)
vcf_all<-list.files("01-inputs-vcf", pattern = "*.vcf", recursive = TRUE, full.names = TRUE)

attackVCF <- function(x) {
      vcfs<-readVcf(x)
      vcfs<-data.frame("chrom_pos"=rownames(data.frame(info(vcfs))))
      return(vcfs)
    }

samples_vcf<-unique(titan_segs$uniq_id)
#strsplit(vcf_all, "/")[[3]][2]

## make position files from vcfs 
for(ii in 1:length (samples_vcf)){

  vcf_per_sample<-vcf_all[grep(samples_vcf[ii],vcf_all)]
  vcfs <- lapply(vcf_per_sample, attackVCF)

    for (i in 1:length(vcfs)){
    vcfs[[i]]<-cbind(vcfs[[i]],vcf_per_sample[i])
    }

     vcfs <- do.call("rbind", vcfs)
     vcfs$chrom_position_id<-gsub("_.*$","",vcfs$chrom_pos)
     vcfs$position<-gsub(".*:","",vcfs$chrom_position_id)
     vcfs$chrom<-gsub(":.*$","",vcfs$chrom_position_id)

     vcfs$specimen<-gsub("01-inputs-vcf/", "",gsub(".chrchr1.vcf","",vcfs[,2]))
     vcfs$sample_id<-gsub("/.*$","",vcfs$specimen)
     vcfs$vcf_sample_id<-gsub("-.*$","",vcfs$specimen)
     vcfs_final<-vcfs[,c("chrom","position","vcf_sample_id")]
     #vcf_meta<-merge(vcf_final,meta, by.x = "sample_id", by.y = "V1")

     write.table(vcfs_final , file = paste("02-inputs-positions/positions_",samples_vcf[ii],".txt",sep = ""), col.names =  TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
 
}

#### run main part to merge segments for each subect ID
titan_segs$type[titan_segs$Corrected_Call=="NEUT"]<-"neutral"
titan_segs$type[titan_segs$Corrected_Call== "AMP"| titan_segs$Corrected_Call== "GAIN" | titan_segs$Corrected_Call== "HLAMP"]<-"gain"
titan_segs$type[titan_segs$Corrected_Call== "HETD"| titan_segs$Corrected_Call== "HOMD"]<-"loss"


out_res_all<-NULL
for(ss in 1:length(patients_rt)){  ### loop through each patient
#for(ss in 1:2){

    focal_patient<-titan_segs[titan_segs$uniq_id2== patients_rt[ss],]
    #focal_patient$type[focal_patient$Corrected_Copy_Number>2]<-"gain"
    #focal_patient$type[focal_patient$Corrected_Copy_Number<2]<-"loss"
    #focal_patient$type[focal_patient$Corrected_Copy_Number==2]<-"neutral"
    
    
    seg_type<-unique(focal_patient$type)

    sample_vcf<-gsub("_.*$","",patients_rt[ss])
    focal_vcf<-read.table(file = paste ("02-inputs-positions/positions_",sample_vcf,".txt", sep = ""), header = TRUE)

        #out_res_all<-NULL
        for(tt in 1:length(seg_type)){  ### for each patient loop through each CNV type and merge overlapping segments

            focal_segtype<-focal_patient[focal_patient$type==seg_type[tt],]
            focal_segtype_gr<-GRanges(IRanges(start = focal_segtype$Start_Position.bp.,end =  focal_segtype$End_Position),seqnames = focal_segtype$Chromosome)
            focal_merged<-data.frame(reduce(focal_segtype_gr))[,c(1:4)]

                for(gg in 1:nrow(focal_merged)){  ### loop through each merged region

                focal_CNVstate<-focal_merged[gg,]
                ## find rows within the range of focal merged in the patient seg file with similar type
                rows <- focal_segtype$Chromosome == focal_CNVstate$seqnames & focal_segtype$Start_Position.bp. >= focal_CNVstate$start & focal_segtype$End_Position.bp. <= focal_CNVstate$end
                data_range<- focal_segtype[rows, ]
                data_range$seg_length<-abs(data_range$Start_Position.bp. - data_range$End_Position.bp.)
                data_range_max<-data_range[which.max(data_range$seg_length),]

                ### find accurate number of markers for the merged regions from vcf file
                vcfs_focal_merged_inrange<-focal_vcf[(focal_vcf$vcf_sample_id%in%data_range$sample_id & focal_vcf$chrom==focal_CNVstate$seqnames & as.numeric(focal_vcf$position) >= focal_CNVstate$start & as.numeric(focal_vcf$position)  <= focal_CNVstate$end),]
                vcfs_focal_merged_inrange$uniq_posID<-paste(vcfs_focal_merged_inrange$chrom,vcfs_focal_merged_inrange$position, sep = "_")
                snps<-unique(vcfs_focal_merged_inrange$uniq_posID)

                #focal_CNVstate$snp_length_mean <-round(mean(data_range$Length.snp.))
                focal_CNVstate$copy_number_mean <-round(mean(data_range$Corrected_Copy_Number))
                focal_CNVstate$corrected_ratio <-(mean(data_range$Corrected_Ratio))
                focal_CNVstate$copy_number_median <-median(data_range$Corrected_Copy_Number)
                focal_CNVstate$sample_id<-unique(focal_patient$uniq_id2)
                focal_CNVstate$longest_segment<-data_range_max$seg_length
                focal_CNVstate$samples_with_fsegment<-length(unique(data_range$sample_id))
                focal_CNVstate$CNV_state<-seg_type[tt]
                focal_CNVstate$SNPs_length<-length(snps)
                focal_CNVstate$samples_with_fsegmentID<-paste("'",unique(data_range$sample_id),"'",collapse=", ",sep="")

                out_res_all<-rbind(focal_CNVstate,out_res_all)  #### final table for gain
        } ## merged row

    } ## segment type

}     ### ss loop (patientRT)

write.table(out_res_all, file = "out_res_all_merged.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
out_res_all<-read.delim(file = "out_res_all_merged.txt", header = TRUE)

################################################
### identical function (when number of gain and loss segments are equal)
ident <- function(...){
    args <- c(...) 
    if( length( args ) > 2L ){
       #  recursively call ident()
       out <- c( identical( args[1] , args[2] ) , ident(args[-1]))
    }else{
        out <- identical( args[1] , args[2] )
    }    
    return( all( out ) )
}


out_res_patient_s<-unique(out_res_all$sample_id)
out_res_patient<-NULL

for(pp in 1:length(out_res_patient)){

  focal_outres_patient<-out_res_all[out_res_all$sample_id==out_res_patient_s[pp],]
  chroms<-unique(focal_outres_patient$seqnames)

  out_res_chrom<-NULL  ### chrom

  for (rr in 1:length(chroms)){

    #for (rr in 2:2){
    
    focal_chrom<-focal_outres_patient[focal_outres_patient$seqnames  == chroms[rr],]
    
    ##  neutral intersections only
    focal_neutral<-focal_chrom[focal_chrom$CNV_state == "neutral",]
    grippy_neutral<-GRanges(IRanges(start = focal_neutral$start,end =  focal_neutral$end),seqnames = focal_neutral$seqnames)
    focal_merged_neutral<-data.frame(reduce(grippy_neutral))[,c(1:4)]
    focal_merged_neutral_gr<-GRanges(IRanges(start = focal_merged_neutral$start,end =  focal_merged_neutral$end),seqnames = focal_merged_neutral$seqnames)

    ## gain and loss intersections only
    focal_chrom_neutraldroped<-focal_chrom[focal_chrom$CNV_state != "neutral",]
    grippy<-GRanges(IRanges(start = focal_chrom_neutraldroped$start,end =  focal_chrom_neutraldroped$end),seqnames = focal_chrom_neutraldroped$seqnames)
    focal_merged2<-data.frame(reduce(grippy))[,c(1:4)]
    focal_merged2_gr<-GRanges(IRanges(start = focal_merged2$start,end =  focal_merged2$end),seqnames = focal_merged2$seqnames)

    count_overlaps<-countOverlaps(focal_merged2_gr,focal_merged_neutral_gr)
    count_overlaps_neutral<-countOverlaps(focal_merged_neutral_gr,focal_merged2_gr)
  
    if (nrow (focal_merged2)>=1){

      out_res_gl<-NULL
      for(kk in 1:nrow(focal_merged2)){   ### loop through each gain/loss 

        focal_xx<-focal_merged2[kk,]

        if (count_overlaps[kk] >=1 & nrow(focal_merged2)>= 1| count_overlaps[kk] ==0 & nrow(focal_merged2)>= 1) {

        qq_overlapped<-focal_chrom_neutraldroped[focal_chrom_neutraldroped$start >=  focal_xx$start & focal_chrom_neutraldroped$end <= focal_xx$end,]

          if (nrow(qq_overlapped)> 1 & (ident(qq_overlapped$samples_with_fsegment)==TRUE)) {

               qq_overlapped_good<-qq_overlapped[which.max(qq_overlapped$corrected_ratio),]  
            } 
        
             if (nrow(qq_overlapped)> 1 & (ident(qq_overlapped$samples_with_fsegment)==FALSE)){

               qq_overlapped_good<-qq_overlapped[which.max(qq_overlapped$samples_with_fsegment),]  
            }  

             if (nrow(qq_overlapped) == 1){

               qq_overlapped_good<-qq_overlapped
            } 

         }

         out_res_gl<-rbind(qq_overlapped_good,out_res_gl)

      } ###kk loop through each gain/loss

    }  ##(if nrow (focal_merged2)>=1

    #if(count_overlaps_neutral ==0){

    neutral_uniq<-count_overlaps_neutral==0
    beck_neutral<-focal_merged_neutral[neutral_uniq,]

      if (nrow(beck_neutral)>0){
      out_res_neutral<-NULL
        for(jj in 1:nrow(beck_neutral)){

        beck_neutral_focal<-beck_neutral[jj,]
        beck_neutral_focal<-focal_neutral[focal_neutral$start >=  beck_neutral_focal$start & focal_neutral$end <= beck_neutral_focal$end,]
        out_res_neutral<-rbind(beck_neutral_focal,out_res_neutral)
        }

      }
    #out_res_chrom<-rbind(gl_neutral,out_res_chrom)

    #}  ## if count_overlaps_neutral ==0
      gl_neutral<-rbind(out_res_gl,out_res_neutral)
      out_res_chrom<-rbind(gl_neutral,out_res_chrom)
      out_res_chrom<-na.omit(out_res_chrom)
      out_res_chrom_good<-out_res_chrom[!duplicated(out_res_chrom),]
      #out_res_patient<-rbind(out_res_chrom_good,out_res_patient)

  }  ## chrom loop
        out_res_patient<-rbind(out_res_chrom_good,out_res_patient)

}  ## patient loop

all_selected<-out_res_patient[,c("sample_id","seqnames" , "start" , "end","SNPs_length","copy_number_median")]
all_selected$Seg.CN<- (log2(all_selected$copy_number_median) -1 )
all_final<-all_selected[,c(1,2,3,4,7,5)]
all_final$seqnames<-gsub("chr","",all_final$seqnames)

all_final_nodup<-all_final[!duplicated(all_final),]
write.table(all_final_nodup, file = "out_res_all_chroms_merged_regions_per_patient_for_GISTIC.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

   


