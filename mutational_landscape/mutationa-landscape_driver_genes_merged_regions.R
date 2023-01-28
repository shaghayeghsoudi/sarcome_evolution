
############################################################
##### process mutations and run dNdScv fpr Sarcoma samples (all samples) #####
rm(list = ls())
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


mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

mutaion_RT$RTtretment<-ifelse(mutaion_RT$RT_status== "preRT", "naive",
                        ifelse(mutaion_RT$RT_status== "noRT", "naive",
                        ifelse(mutaion_RT$RT_status== "postRT", "treatment",
                        "-")))

mutaion_RT<-mutaion_RT[mutaion_RT$Chromosome!="Start",]
#cols <- c("Impact","Coding")
#mutaion_RT[cols] <- lapply(mutaion_RT[cols], as.character)                   
#mutaion_RT$Impact[mutaion_RT$Coding== "splicing"] <- "Splice sites"
#colnames(mutaion_RT)<-c("sampleID","chr","pos","ref","mut","RT_code","Gene","Impact","RTtretment","RT_staus_code")

mutaion_RT$pos_id<-paste(mutaion_RT$Chromosome,mutaion_RT$Start, sep = "_")
mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id,mutaion_RT$RT_status, sep = "_")


samples<-unique(mutaion_RT$identifier)

out_res<-NULL
for (i in 1:length(samples)){

    mutaion_RT_focal<-mutaion_RT[mutaion_RT$identifier==samples[i],]
    mutaion_RT_focal_nodup<-mutaion_RT_focal[!(duplicated(mutaion_RT_focal$pos_id)),]
    out_res<-rbind(mutaion_RT_focal_nodup,out_res)  
}



############################################################
#######Find recurrent mutations by running dNdSCV ##########
############################################################
#mmm<-mm[,c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2" ,"Tumor_Sample_Barcode")]       
#unique(mmm$Variant_Classification)
# [1] Intron                 IGR                    3'UTR                  5'Flank               
# [5] 3'Flank                Missense_Mutation      Nonsense_Mutation      RNA                   
# [9] Silent                 5'UTR                  Nonstop_Mutation       Splice_Region         
#[13] Splice_Site            Frame_Shift_Del        Frame_Shift_Ins        In_Frame_Del          
#[17] In_Frame_Ins           Translation_Start_Site

out_res<-out_res[out_res$Start!="777428",]

### adjust maf file format               
out_res$Variant_Classification[out_res$Coding== "splicing"] <- "Splice_Site"
out_res$Variant_Classification[out_res$Coding== "UTR5"] <- "5'UTR"
out_res$Variant_Classification[out_res$Coding== "UTR3"] <- "3'UTR"
out_res$Variant_Classification[out_res$Impact== "synonymous SNV"] <- "Silent"
out_res$Variant_Classification[out_res$Impact== "stopgain"] <- "Nonsense_Mutation"
out_res$Variant_Classification[out_res$Impact== "nonsynonymous SNV"] <- "Missense_Mutation"
out_res$Variant_Classification[out_res$Impact== "frameshift insertion"] <- "Frame_Shift_Ins"

out_res<-out_res[!(out_res$Impact=="unknown"),]
mutaion_RT_good<-out_res%>%drop_na(Variant_Classification)


#mutaion_RT[cols] <- lapply(mutaion_RT[cols], as.character)

mutaion_RT_good<-mutaion_RT_good[,c("identifier","Chromosome", "Start","End","Ref","Alt", "Gene","RT_status", "RT_code","Variant_Classification")]
colnames(mutaion_RT_good)<-c("Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Hugo_Symbol","RT_status", "RT_code","Variant_Classification")
mutaion_RT_good<-mutaion_RT_good[,c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode","RT_status", "RT_code")]

mutaion_RT_good$Variant_Type<-rep("SNP") #### canbe used as maf file 
mutaion_RT_good$Variant_Type[mutaion_RT_good$Variant_Classification=="Frame_Shift_Ins"]<-"INS"

#### save prepared maf like file ####
write.table(mutaion_RT_good,file = "out_res_samples_updated_regions_merged.maf", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
### Run dnds

mutaion_for_dnds<-mutaion_RT_good[,c("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2")]
colnames(mutaion_for_dnds)<-c("sampleID","chr","pos","ref","mut")
#muts<-c("A" ,"C" , "G" , "T")
#
#mutaion_for_dnds<-mutaion_for_dnds[mutaion_for_dnds$ref%in%muts,]
#mutaion_for_dnds<-mutaion_for_dnds[mutaion_for_dnds$mut%in%muts,]

  dndsout = dndscv(mutaion_for_dnds)
  sel_cv = dndsout$sel_cv
  signif_genes = sel_cv[sel_cv$qallsubs_cv<0.01, c("gene_name","qallsubs_cv")]
  write.table(signif_genes ,file ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_sel_cv_fdr0.01_all_samples.table", ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
   

  globaldnds_con<-(dndsout$globaldnds)
  write.table(globaldnds_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_globaldnds_fdr0.01_all_samples.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
  
  annotmuts_con<-(dndsout$annotmuts)
  write.table(annotmuts_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_annotmuts_fdr0.01_all_samples.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE) 
  
  genemuts_con<-(dndsout$genemuts)
  write.table(genemuts_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_genemuts_fdr0.01_all_sample.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
  
  mle_submodel_con<-(dndsout$mle_submodel)
  write.table(mle_submodel_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_mle_submodel_fdr0.01_all_sample.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
  


######################## MafTools ##############
#################################################
#### analyze the data with maftools #############
rm(list = ls())
library("maftools")
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")
mutaion_RT_good<-read.delim(file = "all_samples_updated_converted_like_multi_reginal_combined.maf", header = TRUE)
#mutaion_RT_good<-mutaion_RT_good[!(mutaion_RT_good$Reference_Allele=="-"),]

sar_laml = read.maf(maf = mutaion_RT_good)

#ellis_BD <- shared_clonal_mutations_fix[order(shared_clonal_mutations_fix$sampleID,shared_clonal_mutations_fix$chr,shared_clonal_mutations_fix$pos),]
#ind = which(diff(ellis_BD $Start_Position)==1)
#ellis_consecutive<-ellis_BD [unique(sort(c(ind,ind+1))),]
#toDelete <- seq(1, nrow(ellis_consecutive), 2)
#toDelete_df<-ellis_consecutive[toDelete ,]

#Shows sample summry.
getSampleSummary(sar_laml)
#Shows gene summary.
getGeneSummary(sar_laml)
#shows clinical data associated with samples
getClinicalData(sar_laml)
#Shows all fields in MAF
getFields(sar_laml)#
#Writes maf summary to an output file with basename laml.
#write.mafSummary(maf = sar_laml, basename = 'laml')

###Plotting MAF summary
plotmafSummary(maf = sar_laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

### draw oncolplot
sig_genes<-read.table(file ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_sel_cv_fdr0.01_all_samples.table", header = TRUE)
sar_aml_genes<-sig_genes[!grepl("RBM25",sig_genes$gene_name),]
sar_aml_genes<-sar_aml_genes[!grepl("MYO3A",sar_aml_genes$gene_name),]
sar_aml_genes<-sar_aml_genes[!grepl("OR10G8",sar_aml_genes$gene_name),]
sar_aml_genes<-sar_aml_genes[!grepl("HERC1",sar_aml_genes$gene_name),]



sar_aml_genes<-sar_aml_genes[1:50,1]
#sar_aml_genes<-sar_aml_genes[-c("RBM25","MYO3A","OR10G8","HERC1")])




oncoplot(maf = sar_laml, genes = sar_aml_genes)

#laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
#plotTiTv(res = laml.titv)


maf_meta<-mutaion_RT_good[,c("Tumor_Sample_Barcode","RT_status","RT_code")]
colnames(maf_meta)<-c("Tumor_Sample_Barcode","RT_status","RT_code")
maf_meta<-maf_meta[!(duplicated(maf_meta$Tumor_Sample_Barcode)),]

sar_laml = read.maf(maf = mutaion_RT_good, clinicalData = maf_meta)

#oncoplot(maf = sar_laml, genes = sar_aml_genes, clinicalFeatures = 'RT_status',sortByAnnotation = TRUE,fontSize = 0.4, top = 2)
oncoplot(maf = sar_laml,  clinicalFeatures = 'RT_status',sortByAnnotation = TRUE,fontSize = 0.5, top = 20)


rainfallPlot(maf = sar_laml, detectChangePoints = TRUE, pointSize = 0.4)





##############
#### from nonmerged




mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]

### adjust maf file format               
mutaion_RT$Variant_Classification[mutaion_RT$Coding== "splicing"] <- "Splice_Site"
mutaion_RT$Variant_Classification[mutaion_RT$Coding== "UTR5"] <- "5'UTR"
mutaion_RT$Variant_Classification[mutaion_RT$Coding== "UTR3"] <- "3'UTR"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "synonymous SNV"] <- "Silent"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "stopgain"] <- "Nonsense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "nonsynonymous SNV"] <- "Missense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "frameshift insertion"] <- "Frame_Shift_Ins"

mutaion_RT<-mutaion_RT[!(mutaion_RT$Impact=="unknown"),]
mutaion_RT_good<-mutaion_RT%>%drop_na(Variant_Classification)



mutaion_RT_good<-mutaion_RT_good[,c("sample_id","Chromosome", "Start","End","Ref","Alt", "Gene","RT_status", "RT_code","Variant_Classification")]
colnames(mutaion_RT_good)<-c("Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Hugo_Symbol","RT_status", "RT_code","Variant_Classification")
mutaion_RT_good<-mutaion_RT_good[,c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode","RT_status", "RT_code")]

mutaion_RT_good$Variant_Type<-rep("SNP") #### canbe used as maf file 
mutaion_RT_good$Variant_Type[mutaion_RT_good$Variant_Classification=="Frame_Shift_Ins"]<-"INS"

#### save prepared maf like file ####
write.table(mutaion_RT_good,file = "out_res_samples_updated_regions_NOTmerged.maf", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
### Run dnds

mutaion_for_dnds<-mutaion_RT_good[,c("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2")]
colnames(mutaion_for_dnds)<-c("sampleID","chr","pos","ref","mut")
#muts<-c("A" ,"C" , "G" , "T")
#


## maftools


mutaion_RT_good<-read.delim(file = "out_res_samples_updated_regions_merged.maf", header = TRUE)
#mutaion_RT_good<-mutaion_RT_good[!(mutaion_RT_good$Reference_Allele=="-"),]
mutaion_RT_good$Variant_Type<-rep("SNP") #### canbe used as maf file 
mutaion_RT_good$Variant_Type[mutaion_RT_good$Variant_Classification=="Frame_Shift_Ins"]<-"INS"

sar_aml_genes<-signif_genes[1:20,1]

oncoplot(maf = sar_laml, genes =   sar_aml_genes)



maf_meta<-mutaion_RT_good[,c("Tumor_Sample_Barcode","RT_status","RT_code")]
colnames(maf_meta)<-c("Tumor_Sample_Barcode","RT_status","RT_code")
maf_meta<-maf_meta[!(duplicated(maf_meta$Tumor_Sample_Barcode)),]

sar_laml = read.maf(maf = mutaion_RT_good, clinicalData = maf_meta)

#oncoplot(maf = sar_laml, genes = sar_aml_genes, clinicalFeatures = 'RT_status',sortByAnnotation = TRUE,fontSize = 0.4, top = 2)
oncoplot(maf = sar_laml,  clinicalFeatures = 'RT_status',sortByAnnotation = TRUE,fontSize = 0.5, top = 20)


rainfallPlot(maf = sar_laml, detectChangePoints = TRUE, pointSize = 0.4)
