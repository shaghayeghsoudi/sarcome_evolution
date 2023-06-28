
### this is new file

### mutational landscape (SNV and CNV-GISTIC)-Oncoplot
### original author: Shaghayegh Soudi
### May 2023
rm(list = ls())
library("tidyverse")
library("plyr")
library("maftools")
library("stringr")
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")

meta$RTstatus<-ifelse(meta$sequenceofsamplevRT== "beforeRT", "noRT",
                        ifelse(meta$sequenceofsamplevRT== "afterRT", "postRT",
                        ifelse(meta$sequenceofsamplevRT== "nopreopRT", "noRT",
                        "-")))


meta_good<-meta%>%mutate(unique_sample_id=gsub("_.*$","",sampleid))%>%
  mutate(identifier=paste(unique_sample_id,RTstatus, sep = "_"))


######################################
### load and read montecarlo files ###
######################################
filenames_ssm_cf <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/MRS_evolution_manuscript_data/cfDNA_analysis/montecarlo/varpossnvs",pattern="*.txt", full.names = TRUE)
attackStats_montecarlo <- lapply(filenames_ssm_cf,function(x) {
     read.delim(x, header=TRUE, sep = "\t",)
     })

### add sample name as a column
for (i in 1:length(attackStats_montecarlo)){
    attackStats_montecarlo[[i]]<-cbind(attackStats_montecarlo[[i]],filenames_ssm_cf[i])
    }
aa_montecarlo<- do.call("rbind", attackStats_montecarlo) 
aa_montecarlo %>% separate (SAMPLE,c("sample","sample_id","general_info"))



montecarlo<-aa_montecarlo %>% separate (SAMPLE,c("sample","sample_id","general_info"))
montecarlo_tumour<-montecarlo[grepl("T",montecarlo$general_info),]   ### subset muttaions only to tumour mutations


montecarlo_tumour$general_info<-gsub("T","",montecarlo_tumour$general_info)
montecarlo_tumour$meta_id<-paste(montecarlo_tumour$sample_id,montecarlo_tumour$general_info, sep = "_")
montecarlo_tumour_meta<-merge(montecarlo_tumour,meta, by.x ="meta_id" , by.y = "sampleid")

montecarlo_tumour_meta$chrom_pos<-paste(montecarlo_tumour_meta$CHR,montecarlo_tumour_meta$POSITION , sep = "-")
montecarlo_tumour_meta$monto_ID<-paste(montecarlo_tumour_meta$identifier, montecarlo_tumour_meta$chrom_pos,  sep = "-")
montecarlo_tumour_meta$monto_ID<-gsub("chr","",montecarlo_tumour_meta$monto_ID)

##match_colB<-c("identifier","CHR", "POSITION")
#montecarlo_tumour_meta$monto_ID <- apply( montecarlo_tumour_meta[ , match_colB] , 1 , paste0 , collapse = "-" )

montecarlo_tumour_meta_dedup<-montecarlo_tumour_meta[!duplicated(montecarlo_tumour_meta$monto_ID),]
monto_df<-montecarlo_tumour_meta_dedup[,c("meta_id","identifier","monto_ID")]  ## final table to be merged with variants raw

###################################
### load and read all maf files ###
###################################
filenames_CN <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/selected",pattern="*.filtered.variants.oxomerge.final.txt_selected", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=FALSE, sep = "\t",)
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
shared_clonal_mutations_fix$sample_id2<-sub('.*/\\s*', '', shared_clonal_mutations_fix$sample_id)
mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id2", by.y = "sampleid")

mutaion_RT$chrom_pos<-paste(mutaion_RT$Chromosome , mutaion_RT$Start , sep = "-")
mutaion_RT$var_id<-paste(mutaion_RT$identifier,mutaion_RT$chrom_pos , sep = "-")

mutaion_RT$Impact[mutaion_RT$Coding== "splicing"] <- "Splice_Site"
mutaion_RT$Impact[mutaion_RT$Coding== "UTR5"] <- "5'UTR"
mutaion_RT$Impact[mutaion_RT$Coding== "UTR3"] <- "3'UTR"


mutaion_RT$Impact <- gsub("synonymous SNV", "Silent",
           gsub("nonsynonymous SNV", "Missense_Mutation",
           gsub("stopgain", "Nonsense_Mutation", 
           gsub("stoploss", "Nonsense_Mutation",
           gsub("frameshift deletion", "Frame_Shift_Del",
           gsub("frameshift insertion", "Frame_Shift_Ins",
           gsub("unknown", "Splice_Site",
           mutaion_RT$Impact)))))))

mutaion_RT$Variant_Type<-rep("SNP")

mutaion_RT_dedup<-mutaion_RT[!duplicated(mutaion_RT$var_id),]    ### remove duplicated rows

##############################
### merge variant files with monto files ###
##############################
mutaion_RT_monto<-merge(mutaion_RT_dedup,monto_df,by.x="var_id",by.y="monto_ID")

mutaion_RT_monto_good<-mutaion_RT_monto[,c("Gene","Chromosome","Start","End","Impact","Ref", "Alt","identifier.x","RTstatus","Variant_Type")]
mutaion_RT_monto_good<-mutaion_RT_monto_good[mutaion_RT_monto_good$Start!="Start",]

colnames(mutaion_RT_monto_good)<-c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode","RT_status","Variant_Type")

bad_impact<-c("Coding_info")
mutaion_RT_monto_good<-mutaion_RT_monto_good[!(mutaion_RT_monto_good$Variant_Classification%in%bad_impact),]
write.table(mutaion_RT_monto_good,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/pre_postRT_samples_updated_converted_like_updatedmonto.maf",col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)



sar_aml_genes<-c("TP53","MUC16","TNC","CLTC","NBEA","ATRX","RBM10","COL2A1")


#cols <- c("Impact","Coding")
#mutaion_RT[cols] <- lapply(mutaion_RT[cols], as.character)                   
#mutaion_RT$Impact[mutaion_RT$Coding== "splicing"] <- "Splice sites"
#colnames(mutaion_RT)<-c("sampleID","chr","pos","ref","mut","RT_code","Gene","Impact","RTtretment","RT_staus_code")

##################################
#######################################
###### plot variants with maf tools ###
colnames(meta)[9]<-"Tumor_Sample_Barcode"
mutaion_RT_good_monto<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/pre_postRT_samples_updated_converted_like_updatedmonto.maf", header = TRUE)  ### I changed location to GISTIC analysis folder
mutaion_RT_good$Variant_Classification<-gsub("Splice sites","Splice_Site",mutaion_RT_good$Variant_Classification)


vc_cols = RColorBrewer::brewer.pal(n = 6, name = 'Accent')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'Splice_Site'
)

#Color coding for FAB classification; try getAnnotations(x = laml) to see available annotations.
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6")

#vc_nonsyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")
vc_nonsyn = c( "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation", "Missense_Mutation")
vc_nonsyn = c(vc_nonsyn, "Silent")



ordered<-c("SRC125_noRT","SRC127_noRT","SRC130_noRT","SRC167_noRT","SRC168_noRT","SRC169_noRT","SRC170_noRT","SRC171_noRT","SRC172_noRT","SRC173_noRT","TB13092_noRT","TB13712_noRT","TB13959_noRT","TB8016_noRT","TB9051_noRT","TB9573_noRT","SRC125_postRT","SRC127_postRT","SRC150_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB11985_postRT","TB12052_postRT","TB22446_postRT","TB9051_postRT")


laml = read.maf(maf = mutaion_RT_good_monto,
                clinicalData = meta,
                verbose = FALSE,vc_nonSyn = vc_nonsyn)

#Shows sample summry.
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#Shows all fields in MAF
getFields(laml) 
#Writes maf summary to an output file with basename laml.
#write.mafSummary(maf = sar_laml, basename = 'laml')


oncoplot(maf = laml, draw_titv = TRUE ,
        colors = vc_cols,genes = sar_aml_genes,
        clinicalFeatures =c("RTstatus","Gender"),
        sortByAnnotation = TRUE,drawRowBar = TRUE,
        removeNonMutated= FALSE,barcode_mar = 10,
        sortByMutation=TRUE,annoBorderCol="white",
        fontSize=0.5,bgCol="white",
        #showTumorSampleBarcodes = TRUE,
        #SampleNamefontSize = 1,
        sampleOrder=ordered)

##################################################
### compare two cohorts (before v. after radiation)
noRT<-mutaion_RT_good[grepl("noRT",mutaion_RT_good$RT_status),]
postRT<-mutaion_RT_good[grepl("postRT",mutaion_RT_good$RT_status),]

noRT.lam = read.maf(maf = noRT)
postRT.lam = read.maf(maf = postRT)

#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = noRT.lam, m2 = postRT.lam , m1Name = 'noRT', m2Name = 'postRT', minMut = 7)
print(pt.vs.rt)
#########################
####### plot GISTIC #####

#########
#par(mar = c(10, 5, 5, 5))
all_lesions<-read.delim("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/all_lesions.conf_75.txt")

### select colnames with sample iD
gistic_samples<-all_lesions %>%
  select(matches('SRC|TB')) %>%
  colnames()

meta_gistic<-meta[meta$sampleid%in%gistic_samples,]
meta_gistic<-meta_gistic[,"Tumor_Sample_Barcode"]


colnames(all_lesions)[10:36]<-meta_gistic
all_lesions <- all_lesions [1: ncol(all_lesions)-1 ]
write.table(all_lesions,file = "~/Desktop/out_gistic/all_legions_updated.txt", col.names = TRUE, row.nam = FALSE, sep = "\t", quote = FALSE)



vc_cols2 = c("#0484e6","tomato")
names(vc_cols2) = c(
  'Del',
  'Amp'
)


laml.gistic_del <- readGistic(gisticAllLesionsFile ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/all_lesions.conf_75.txt", 
                                                  #gisticAmpGenesFile = "~/Desktop/out_gistic/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/amp_genes.conf_75.txt", 
                                                  gisticDelGenesFile = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/del_genes.conf_75.txt", 
                                                  cnLevel = "all", gisticScoresFile = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/scores.gistic")



gisticOncoPlot(gistic = laml.gistic_del,  
    removeNonAltered = FALSE,
    #clinicalData = getClinicalData(x = laml),
    #clinicalFeatures = "RTstatus",
    sortByAnnotation = TRUE,
    colors = vc_cols2,bgCol="white",
    fontSize=0.5,
    showTumorSampleBarcodes = TRUE, 
    SampleNamefontSize = 0.8,
    sampleOrder=ordered)


##### Amplifications (temporray path, works)
laml.gistic <- readGistic(gisticAllLesionsFile ="~/Desktop/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/all_lesions.conf_75.txt", 
                          gisticAmpGenesFile = "~/Desktop/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/amp_genes.conf_75.txt", 
                          #gisticDelGenesFile = "~/Desktop/out_gistic/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/del_genes.conf_75.txt", 
                          cnLevel = "all", gisticScoresFile = "~/Desktop/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/scores.gistic")

bands_amp<-c("AP_2:1p32.1","AP_3:1p31.2","AP_5:1q21.1","AP_6:3p12.1","AP_7:5p15.33", "AP_8:7p22.1","AP_9:17p11.2" ,"AP_10:20q11.22", "AP_12:20q13.2","AP_13:20q13.33", "AP_15:22q12.3")
gisticOncoPlot(gistic = laml.gistic,  removeNonAltered = FALSE,clinicalData = getClinicalData(x = laml),
    #clinicalFeatures = "RTstatus",
    sortByAnnotation = TRUE,
    colors = vc_cols2,fontSize=0.5,
    bgCol="white",bands=bands_amp,
    showTumorSampleBarcodes = TRUE, 
    #SampleNamefontSize = 0.6,
    sampleOrder=ordered)


