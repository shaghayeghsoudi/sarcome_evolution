

### make MUTSIG2CV input file
############################################################
##### process mutations and run dNdScv fpr Sarcoma samples (all samples) #####
rm(list = ls())
library("tidyverse")
library("plyr")

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")

## all samples
meta<-read.table(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)


### load and read all maf files
filenames_CN <- list.files("filtered_variants/selected",pattern="*.filtered.variants.oxomerge.final.txt_selected", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=FALSE, sep = "\t")
     })

### add sample name as a column
for (i in 1:length(attackStats_CN)){
    attackStats_CN[[i]]<-cbind(attackStats_CN[[i]],filenames_CN[i])
    }
aa <- do.call("rbind", attackStats_CN) 

aa$sample_id<-gsub("filtered_variants/selected/", "",gsub(".filtered.variants.oxomerge.final.txt_selected","",aa[,10]))
shared_clonal_mutations_fix<-aa[-1,c(-10,-8)]
colnames(shared_clonal_mutations_fix)<-c("Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","sample_id")


#shared_clonal_mutations_fix<-shared_clonal_mutations%>%drop_na(sample_id)
shared_clonal_mutations_fix$Chromosome<-gsub("chr","",shared_clonal_mutations_fix$Chromosome)
shared_clonal_mutations_fix$pos_id<-paste(shared_clonal_mutations_fix$Chromosome,shared_clonal_mutations_fix$Start, sep = "_")

## all samples
meta<-read.table(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)



mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")
mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]

mut_777<-mutaion_RT[mutaion_RT$Start=="777428",]
mut_777<-mut_777[!duplicated(mut_777$pos_id),]

mutaion_RT<-rbind(mutaion_RT,mut_777)

mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

mutaion_RT$RTtretment<-ifelse(mutaion_RT$RT_status== "preRT", "naive",
                        ifelse(mutaion_RT$RT_status== "noRT", "naive",
                        ifelse(mutaion_RT$RT_status== "postRT", "treatment",
                        "-")))

mutaion_RT<-mutaion_RT[mutaion_RT$Chromosome!="Start",]


mutaion_RT$pos_id<-paste(mutaion_RT$Chromosome,mutaion_RT$Start, sep = "_")
mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id,mutaion_RT$RT_status, sep = "_")


samples<-unique(mutaion_RT$identifier)

out_res<-NULL
for (i in 1:length(samples)){

    mutaion_RT_focal<-mutaion_RT[mutaion_RT$identifier==samples[i],]
    mutaion_RT_focal_nodup<-mutaion_RT_focal[!(duplicated(mutaion_RT_focal$pos_id)),]
    out_res<-rbind(mutaion_RT_focal_nodup,out_res)  
}

out_res$Variant_Classification[out_res$Impact== "synonymous SNV"] <- "silent"
out_res$Variant_Classification[out_res$Impact== "nonsynonymous SNV"] <- "Missense_Mutation"
out_res$Variant_Classification[out_res$Impact== "frameshift insertion"] <- "Frame_Shift_Ins"
out_res$Variant_Classification[out_res$Impact== "frameshift deletion"] <- "Frame_Shift_Del"
out_res$Variant_Classification[out_res$Impact== "stopgain"] <- "Nonsense_Mutation"
out_res$Variant_Classification[out_res$Impact== "stoploss"] <- "Nonsense_Mutation"
out_res$Variant_Classification[out_res$Impact== "nonframeshift insertion"] <- "Non_Frame_Shift_Ins"
out_res$Variant_Classification[out_res$Impact== "nonframeshift deletion"] <- "Non_Frame_Shift_Del"
out_res<-out_res[!(out_res$Impact=="unknown"),]
out_res<-out_res[!(out_res$Impact=="."),]



out_res$classification[out_res$Impact== "synonymous SNV"] <- "SNP"
out_res$classification[out_res$Impact== "nonsynonymous SNV"] <- "SNP"
out_res$classification[out_res$Impact== "frameshift insertion"] <- "INS"
out_res$classification[out_res$Impact== "frameshift deletion"] <- "DEL"
out_res$classification[out_res$Impact== "stopgain"] <- "SNP"
out_res$classification[out_res$Impact== "stoploss"] <- "SNP"
out_res$classification[out_res$Impact== "nonframeshift insertion"] <- "INS"
out_res$classification[out_res$Impact== "nonframeshift deletion"] <- "DEL"


out_res_final<-out_res[,c("Chromosome","Start","Gene","identifier","Ref","Alt","Variant_Classification","classification")]
colnames(out_res_final)<-c("chr","pos","gene","patient","ref_allele","newbase","type","classification")
write.table(out_res_final, file = "out_res_MutSig2CV_input_merged_regions.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)