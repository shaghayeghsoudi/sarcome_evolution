

#### make inputs for MutationalSignature from deep sequencing and run the analysis #####
#### muttaional signature analysis for "deep sequencing" data with merged regions ######

## deep sequencing data 
rm(list = ls())

library(stringr)
library(dplyr)
library(reshape2)
library(kableExtra)
library(ggplot2)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(pheatmap)

# Load mutSignatures
library(mutSignatures)
# prep hg19
hg19 <- BSgenome.Hsapiens.UCSC.hg19
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")

## load meta data (available samples)
meta<-read.table(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)


### make input file for mutsignatures ###
### load and read all maf files
filenames_CN <- list.files("cfDNA_analysis/snvs",pattern="*_cfdna.All.varpos.snvs.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

### add sample name as a column
for (i in 1:length(attackStats_CN)){
    attackStats_CN[[i]]<-cbind(attackStats_CN[[i]],filenames_CN[i])
    }
aa <- do.call("rbind", attackStats_CN) 
aa$sample_id<-str_split_fixed(aa$SAMPLE, '_', 3)[,2]

aa$CHR<-gsub("chr","",aa$CHR)
aa$pos_id<-paste(aa$CHR,aa$POSITION, sep = "_")



cf_mutations<-aa[,-c(1,20)]
cf_mutations$sample_id<-gsub("-","_",cf_mutations$sample_id)
cf_mutations_fix<-cf_mutations%>%separate(sample_id,c("uniq_identifier","cf_status"))
cf_mutations_tumour<-cf_mutations_fix[grepl("T",cf_mutations_fix$cf_status),]
cf_mutations_tumour$cf_status<-gsub("T","",cf_mutations_tumour$cf_status)
cf_mutations_tumour$sample_id<-paste(cf_mutations_tumour$uniq_identifier,cf_mutations_tumour$cf_status, sep = "_")

cf_mutations_tumour_good_samples<-merge(cf_mutations_tumour,meta, by.x = "sample_id", by.y = "sample_id")


#colnames(shared_clonal_mutations_fix)<-c("Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","sample_id")

cf_mutations_tumour_good_samples$RT_status<-gsub("noRT-1","noRT",cf_mutations_tumour_good_samples$RT_status)
cf_mutations_tumour_good_samples$RT_status<-gsub("noRT-2","noRT",cf_mutations_tumour_good_samples$RT_status)

cf_mutations_tumour_good_samples$RTtretment<-ifelse(cf_mutations_tumour_good_samples$RT_status== "preRT", "naive",
                        ifelse(cf_mutations_tumour_good_samples$RT_status== "noRT", "naive",
                        ifelse(cf_mutations_tumour_good_samples$RT_status== "postRT", "treatment",
                        "-")))


cf_mutations_tumour_good_samples$identifier<-paste(cf_mutations_tumour_good_samples$unique_sample_id,cf_mutations_tumour_good_samples$RT_status, sep = "_")

sample_regions<-unique(cf_mutations_tumour_good_samples$sample_id)

out_res_cf<-NULL
for(kk in 1:length(sample_regions)){

    cf_muttaions_focal<-cf_mutations_tumour_good_samples[cf_mutations_tumour_good_samples$sample_id==sample_regions[kk],]
    cf_muttaions_focal_TP<-cf_muttaions_focal[cf_muttaions_focal$TUMOR_PERCENT>=5,]
    out_res_cf<-rbind(cf_muttaions_focal_TP,out_res_cf)
}



samples<-unique(out_res_cf$identifier)

out_res<-NULL
for (i in 1:length(samples)){

    mutaion_RT_focal<-out_res_cf[out_res_cf$identifier==samples[i],]
    mutaion_RT_focal_nodup<-mutaion_RT_focal[!(duplicated(mutaion_RT_focal$pos_id)),]
    out_res<-rbind(mutaion_RT_focal_nodup,out_res)  
}


write.table(out_res, file = "outres_cfdna.All.varpos.snvs.merged_regions_all_samples_threshold_5tumourpercent.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

out_res_good<-out_res[,c("CHR", "POSITION","NORMAL_ALLELE","TUMOR_ALLELE","identifier")]
colnames(out_res_good)<-c("CHROM","POS","REF","ALT","SAMPLEID")

### drop samples with very low muttaional burden
count<-table(out_res_good$SAMPLEID)
out_res_good_samples<-out_res_good[out_res_good$SAMPLEID %in% names(count[count > 20]), ]


out_res_good_samples$CHROM<-paste("chr",out_res_good_samples$CHROM,sep = "")

## write input file for mutsignatures ###
write.table(out_res_good_samples, file = "outres_input_mutsignatures_merged_regions_all_cfDNA_samples_threshold_5tumourpercent_tab.table", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


input_sigmutation<-read.table(file = "MutSignatures/outres_input_mutsignatures_merged_regions_all_cfDNA_samples_threshold_5tumourpercent_tab.table", header = TRUE)
### filter out any muttaion that is not a SNV

muts<-c("A","T","C","G")
input_sigmutation_snv<-input_sigmutation[input_sigmutation$REF%in%muts,]
input_sigmutation_snv<-input_sigmutation_snv[input_sigmutation_snv$ALT%in%muts,]
#write.table(input_sigmutation_snv, file = "outres_input_mutsignatures_merged_regions_all_samples_tab.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

########################################################################
########################################################################
########################################################################
#### start with Mutsignatures ####
#### extrcat preRT samples
#y_pre<-input_sigmutation_snv[grepl("preRT",input_sigmutation_snv$SAMPLEID),]
y_pre<-input_sigmutation_snv
#De novo extraction of Mutational Signatures from BLCA samples
# Attach context
y_pre <- attachContext(mutData = y_pre,
                   chr_colName = "CHROM",
                   start_colName = "POS",
                   end_colName = "POS",
                   nucl_contextN = 3,
                   BSGenomeDb = hg19)



# Remove mismatches
y_pre <- removeMismatchMut(mutData = y_pre,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")    
# Compute mutType
y_pre <- attachMutType(mutData = y_pre,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context") 



write.table(y_pre, file = "MutSignatures/outres_compute_mutType_MutSignatures_all_samples_cfDNA_samples_merged_regions.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

######################################################
### look at corrections with particular muttaion type
#age<-read.delim(file= "~/Desktop/age.txt", header = FALSE)
#colnames(age)<-c("sample_id","non","necrosis","age")
#my_y_pre<-y_pre%>%separate(SAMPLEID,c("uniq_sample_id","RT"))
#age_good<-age[age$sample_id%in%my_y_pre$uniq_sample_id,]

#pre_with_age<-merge(my_y_pre,age_good, by.x = "uniq_sample_id", by.y = "sample_id")
#pre_with_age_clock_mut<-pre_with_age[grepl("C>T",pre_with_age$mutType),]
#freq<-data.frame(table(pre_with_age$uniq_sample_id))

#final<-merge(freq, age_good, by.x = "Var1", by.y = "sample_id")

#ggplot(final, aes(x=age, y=Freq)) + 
#geom_point()+
#  geom_smooth(method=lm)

### with ggscatter 
# ggscatter(final, x="age", y="Freq",
#          add = "reg.line",                                 # Add regression line
#          conf.int = TRUE,                                  # Add confidence interval
#         add.params = list(color = "blue",
#                            fill = "lightgray")
#          )+
#  stat_cor(method = "pearson", label.x = 3, label.y = 30)  # Add correlation coefficient

######################################################
######################################################

pre.counts <- countMutTypes(mutTable = y_pre,
                             mutType_colName = "mutType",
                             sample_colName = "SAMPLEID")
# Mutation Counts
print(pre.counts)                  
mouCancer.assess <- prelimProcessAssess(input = pre.counts, approach = "counts")


#sarcoma.counts.data <- getCounts(pre.counts)
x <- getCounts(pre.counts)
#ocd <- as.mutation.counts(sarcoma.counts.data)
xx <- as.mutation.counts(x)
# how many signatures should we extract? 

mouCancer.assess <- prelimProcessAssess(input = xx, approach = "counts")

num.sign <- 4

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
pre.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 500,               # bootstrapping: usually 500-1000
    num_parallelCores = 6)                # total num of cores to use (parallelization)


# Extract new signatures (de-novo mutations)- may take a while
pre.analysis <- 
  decipherMutationalProcesses(input = pre.counts,
                              params = pre.params)


#Downstream analyses and visualization
# Retrieve signatures (results)
pre.sig <- pre.analysis$Results$signatures

# Retrieve exposures (results)
pre.exp <- pre.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)

## => fix the code, loop does not work
for (pp in 1:length(num.sign)){

  pdf(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/denovo_MutSignature_sign_",pp,"_PreRT_merged_regions_boot500_cgDNA.pdf"))
  msigPlot(pre.sig, signature = pp, ylim = c(0, 0.10))
  dev.off()
}


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/barplot_denovo_MutSignature_counts_per_sample_PreRT_merged_regions_boot500_cfDNA.pdf")
msigPlot(pre.exp) + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))
dev.off()
xprt <- coerceObj(x = pre.exp, to = "data.frame") 
#xprt <- coerceObj(x = pre.exp, to = "data.frame") 

write.table(xprt, file = "MutSignatures/cf/out_put_mutsignatures_xprt_denove_signature_count_PreRT_cfDNA.txt")

########################################################
##### costom visualization of denovo exposures #########

denovo_exp<-data.frame(t(xprt))

our_res_exp<-NULL
for (jj in 1:ncol(denovo_exp)){

    denovo_exp_focal<-data.frame(denovo_exp[,jj])
    denovo_exp_focal$sample_id<-rownames(denovo_exp)
    denovo_exp_focal$sig<-rep(colnames(denovo_exp)[jj])
    colnames(denovo_exp_focal)<-c("count","sample_id","sig")
    our_res_exp<-rbind(denovo_exp_focal,our_res_exp)
}


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/barplot_cexp_counts_cfDNA_preRT_sample_annotated.pdf", width = 5, height = 5)
plota<-ggplot(our_res_exp, aes(x = sample_id, y= count,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   


#### heatmap of denovo exposure frequencies ####
freq<-data.frame(table(y_pre$SAMPLEID))
#freq<-freq%>%separate(Var1,c("uniq_sample_id","RT"))


samples<-colnames(xprt)

out_res_xprt<-NULL
for(k in 1:length(samples)){

    xprt_focal<-xprt[,colnames(xprt)==samples[k]]
    focal_count<-freq[freq$Var1==samples[k],]
    table_count<-data.frame("freq"=xprt_focal/focal_count$Freq)
    colnames(table_count)<-samples[k]

    table<-t(table_count)
    #rownames(table)<-c("Sign.01","Sign.02","Sign.03","Sign.04")
    out_res_xprt<-rbind(table,out_res_xprt)

}

out_res_final_xprt<-t(out_res_xprt)
rownames(out_res_final_xprt)<-c("Sign.01","Sign.02","Sign.03","Sign.04")
out_res_final_xprt_df<-data.frame(out_res_final_xprt)

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_clustered.pdf", width = 5, height = 5)
pheatmap(out_res_final_xprt_df,cluster_rows =FALSE)
dev.off()

ordered<-c("SRC125_preRT","SRC127_preRT","SRC167_preRT","SRC169_preRT","SRC170_preRT","SRC171_preRT","TB9051_preRT","SRC125_postRT","SRC127_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB22446_postRT","TB9051_postRT","RC130_noRT","SRC168_noRT" ,"SRC172_noRT","SRC173_noRT", "TB13092_noRT","TB13712_noRT","TB13959_noRT","TB9573_noRT")
out_res_final_xprt_ordered<-out_res_final_xprt_df[order(match(colnames(out_res_final_xprt_df), ordered))]

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_ordered_Unclustered.pdf", width = 5, height = 5)
pheatmap(out_res_final_xprt_ordered,cluster_cols=FALSE)
dev.off()


##### compare if there is any difference between signatures among different pre,post and noRT samples ###
sig1_deno<-data.frame("proportion"=out_res_final_xprt[1,])
sig1_deno$sample_id<-rownames(sig1_deno)
sig1_deno<-sig1_deno%>%separate(sample_id,c("uniq_sample_id","RT"))


sig2_deno<-data.frame("proportion"=out_res_final_xprt[2,])
sig2_deno$sample_id<-rownames(sig2_deno)
sig2_deno<-sig2_deno%>%separate(sample_id,c("uniq_sample_id","RT"))


sig3_deno<-data.frame("proportion"=out_res_final_xprt[3,])
sig3_deno$sample_id<-rownames(sig3_deno)
sig3_deno<-sig3_deno%>%separate(sample_id,c("uniq_sample_id","RT"))


sig4_deno<-data.frame("proportion"=out_res_final_xprt[4,])
sig4_deno$sample_id<-rownames(sig4_deno)
sig4_deno<-sig4_deno%>%separate(sample_id,c("uniq_sample_id","RT"))

### box plots of difference betwen groups in each signature ###

plot_sig1<-ggplot(sig1_deno, aes(x=RT, y=proportion, fill = RT)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="denovo_exposure_sign01")


plot_sig2<-ggplot(sig2_deno, aes(x=RT, y=proportion, fill = RT)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
                labs(title="denovo_exposure_sign02")


plot_sig3<-ggplot(sig3_deno, aes(x=RT, y=proportion, fill = RT)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
                labs(title="denovo_exposure_sign03")


plot_sig4<-ggplot(sig4_deno, aes(x=RT, y=proportion, fill = RT)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
                labs(title="denovo_exposure_sign04")


library("gridExtra")

grid.arrange(arrangeGrob(plot_sig1, plot_sig2, ncol = 2),                             # First row with one plot spaning over 2 columns
             arrangeGrob(plot_sig3, plot_sig4, ncol = 2), # Second row with 2 plots in 2 different columns
             nrow = 2)                       # Number of rows
###############################################################
###############################################################
# Retrieve COSMIC signatures from online repo, and then subset
###############################################################
###############################################################
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
cosmx<-cosmix[c(1, 2, 5,13)]

# match OVcar and COSMIC signatures
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_similarities_denovo_cosmic_cfDNA_preRT.pdf", width = 5, height = 5)
mSign.sar <- matchSignatures(mutSign = pre.sig, reference = cosmix)
print(mSign.sar$plot)
dev.off()


blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = pre.sig)

blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.1x_df <- coerceObj(x = blca.exp.1x, to = "data.frame") 
write.table(blca.exp.1x_df, file = "MutSignatures/cf/out_put_mutsignatures_blca.exp.1x_cosmic_signature_count_preRT_cfDNA.txt")


blca.exp.2x <- blca.expo2$Results$count.result
blca.exp.2x_df <- coerceObj(x = blca.exp.2x, to = "data.frame") 
write.table(blca.exp.2x_df, file = "MutSignatures/cf/out_put_mutsignatures_blca.exp.2x_denove_signature_count_preRT_cfDNA.txt")


# Plot exposures
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/barplot_cosmic_signatures_bp1_counts_cfDNA_preRTsample.pdf")
bp1 <- msigPlot(blca.exp.1x, main = "BLCA | COSMIC sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))
print(bp1)
dev.off()

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/barplot_bp2_signatures_counts_cfDNA_preRTsample.pdf")
bp2 <- msigPlot(blca.exp.2x, main = "BLCA | pre-blca sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))
print(bp2)
dev.off()


# Visualize
grid.arrange(bp1, bp2, ncol = 2)


#### custom plots of cosmic signatures
cosmic_exp<-data.frame(t(blca.exp.1x_df))

out_res_cos<-NULL
for (jj in 1:ncol(cosmic_exp)){

    cosmic_exp_focal<-data.frame(cosmic_exp[,jj])
    cosmic_exp_focal$sample_id<-rownames(cosmic_exp)
    cosmic_exp_focal$sig<-rep(colnames(cosmic_exp)[jj])
    colnames(cosmic_exp_focal)<-c("count","sample_id","sig")
    out_res_cos<-rbind(cosmic_exp_focal,out_res_cos)
}


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/barplot_cosmic_signatures_counts_cfDNA_sample_annotated.pdf", width = 5, height = 5)
plotb<-ggplot(out_res_cos, aes(x = sample_id, y= count,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plotb)
dev.off()   



#### heatmap of denovo exposure frequencies ####
freq<-data.frame(table(y_pre$SAMPLEID))
#freq<-freq%>%separate(Var1,c("uniq_sample_id","RT"))


samples<-colnames(blca.exp.1x_df)

out_res_xprtcos<-NULL
for(k in 1:length(samples)){

    xprtcos_focal<-blca.exp.1x_df[,colnames(blca.exp.1x_df)==samples[k]]
    focal_count<-freq[freq$Var1==samples[k],]
    table_count<-data.frame("freq"=xprtcos_focal/focal_count$Freq)
    colnames(table_count)<-samples[k]

    table<-t(table_count)
    #rownames(table)<-c("Sign.01","Sign.02","Sign.03","Sign.04")
    out_res_xprtcos<-rbind(table,out_res_xprtcos)

}

out_res_final_xprtcos<-t(out_res_xprtcos)
rownames(out_res_final_xprtcos)<-c("Cosmic.01","Cosmic.02","Cosmic.03","Cosmic13")
out_res_final_xprtcos_df<-data.frame(out_res_final_xprtcos)

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_Cosmic_signatures_frequencies_cfDNA_sample_annotated_clustered.pdf", width = 5, height = 5)
pheatmap(out_res_final_xprtcos_df,cluster_rows =FALSE)
dev.off()

ordered<-c("SRC125_preRT","SRC127_preRT","SRC167_preRT","SRC169_preRT","SRC170_preRT","SRC171_preRT","TB9051_preRT","SRC125_postRT","SRC127_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB22446_postRT","TB9051_postRT","RC130_noRT","SRC168_noRT" ,"SRC172_noRT","SRC173_noRT", "TB13092_noRT","TB13712_noRT","TB13959_noRT","TB9573_noRT")
out_res_final_xprt_ordered<-out_res_final_xprt_df[order(match(colnames(out_res_final_xprt_df), ordered))]

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_ordered_Unclustered.pdf", width = 5, height = 5)
pheatmap(out_res_final_xprt_ordered,cluster_cols=FALSE)
dev.off()


################################################################
################################################################
#### look at the correlation between signatures and age ########
##### => Needs to be optimized <= ##############################
freq<-data.frame(table(y_pre$SAMPLEID))
freq<-freq%>%separate(Var1,c("uniq_sample_id","RT"))

#age_good
colnames(age_good)[1]<-"uniq_sample_id"
freq_age<-merge(freq,age_good, by="uniq_sample_id", all.x=T)

### load signatures ###
#write.table(xprt, file = "MutSignatures/out_put_mutsignatures_xprt_denove_signature_count_All_RT.txt")
xprt<-read.table(file = "MutSignatures/out_put_mutsignatures_xprt_denove_signature_count_All_RT.txt", header = TRUE)


vv<-data.frame(t(xprt))

aa<-data.frame(vv[,1])
aa$sample_id<-rownames(vv)
aa$sig<-rep("sig.1")
colnames(aa)<-c("count","sample_id","sig")


bb<-data.frame(vv[,2])
bb$sample_id<-rownames(vv)
bb$sig<-rep("sig.2")
colnames(bb)<-c("count","sample_id","sig")


cc<-data.frame(vv[,3])
cc$sample_id<-rownames(vv)
cc$sig<-rep("sig.3")
colnames(cc)<-c("count","sample_id","sig")



dd<-data.frame(vv[,4])
dd$sample_id<-rownames(vv)
dd$sig<-rep("sig.4")
colnames(dd)<-c("count","sample_id","sig")

all_me<-rbind(aa,bb,cc,dd)


####
pdf(file = "~/Desktop/scatter_cor_denovo-exposure_sig1_preRT.pdf")
all_me_sig1<-all_me[all_me$sig=="sig.1",]
all_me_sig1<-all_me_sig1%>%separate(sample_id,c("uniq_sample_id","RT"))
all_me_sig1_final<-merge(all_me_sig1,age_good, by="uniq_sample_id", all.x=T)
all_me_sig1_final_pre<-all_me_sig1_final[all_me_sig1_final$RT=="preRT",]
plotaa<-ggscatter(all_me_sig1_final_pre, x="age", y="count",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
         add.params = list(color = "red",
                            fill = "lightgrey")
          )+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)  # Add correlation coefficient
print(plotaa)
dev.off()


pdf(file = "~/Desktop/scatter_cor_denovo-exposure_sig2_preRT.pdf")
all_me_sig2<-all_me[all_me$sig=="sig.2",]
all_me_sig2<-all_me_sig2%>%separate(sample_id,c("uniq_sample_id","RT"))
all_me_sig2_final<-merge(all_me_sig2,age_good, by="uniq_sample_id", all.x=T)
all_me_sig2_final_pre<-all_me_sig2_final[all_me_sig2_final$RT=="preRT",]
plotbb<-ggscatter(all_me_sig2_final_pre, x="age", y="count",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
         add.params = list(color = "red",
                            fill = "lightgrey")
          )+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)  # Add correlation coefficient
print(plotbb)
dev.off()



pdf(file = "~/Desktop/scatter_cor_denovo-exposure_sig3_preRT.pdf")
all_me_sig3<-all_me[all_me$sig=="sig.3",]
all_me_sig3<-all_me_sig3%>%separate(sample_id,c("uniq_sample_id","RT"))
all_me_sig3_final<-merge(all_me_sig3,age_good, by="uniq_sample_id", all.x=T)
all_me_sig3_final_pre<-all_me_sig3_final[all_me_sig3_final$RT=="preRT",]
plotaa<-ggscatter(all_me_sig3_final_pre, x="age", y="count",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
         add.params = list(color = "red",
                            fill = "lightgrey")
          )+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)  # Add correlation coefficient
print(plotaa)
dev.off()


pdf(file = "~/Desktop/scatter_cor_denovo-exposure_sig4_preRT.pdf")
all_me_sig4<-all_me[all_me$sig=="sig.4",]
all_me_sig4<-all_me_sig4%>%separate(sample_id,c("uniq_sample_id","RT"))
all_me_sig4_final<-merge(all_me_sig4,age_good, by="uniq_sample_id", all.x=T)
all_me_sig4_final_pre<-all_me_sig4_final[all_me_sig4_final$RT=="preRT",]
plotaa<-ggscatter(all_me_sig4_final_pre, x="age", y="count",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
         add.params = list(color = "red",
                            fill = "lightgrey")
          )+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)  # Add correlation coefficient
print(plotaa)
dev.off()



##################################################################
##################################################################
#################### postRT ######################################
##################################################################
##################################################################


post<-mutaion_RT[mutaion_RT$RT_status=="postRT",]
y_post<-post[,c("Chromosome","Start","Ref","Alt","sample_id")]
colnames(y_post)<-c("CHROM","POS","REF","ALT","SAMPLEID")


#De novo extraction of Mutational Signatures from BLCA samples
# Attach context
y_post <- attachContext(mutData = y_post,
                   chr_colName = "CHROM",
                   start_colName = "POS",
                   end_colName = "POS",
                   nucl_contextN = 3,
                   BSGenomeDb = hg19)



# Remove mismatches
y_post <- removeMismatchMut(mutData = y_post,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")    
# Compute mutType
y_post <- attachMutType(mutData = y_post,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context") 


#
post.counts <- countMutTypes(mutTable = y_post,
                             mutType_colName = "mutType",
                             sample_colName = "SAMPLEID")
# Mutation Counts
print(post.counts)                  

 mouCancer.assess <- prelimProcessAssess(input = post.counts, approach = "counts")


#sarcoma.counts.data <- getCounts(pre.counts)
x <- getCounts(post.counts)
#ocd <- as.mutation.counts(sarcoma.counts.data)
xx <- as.mutation.counts(x)
# how many signatures should we extract? 

mouCancer.assess <- prelimProcessAssess(input = xx, approach = "counts")

num.sign <- 4

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
post.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 100,               # bootstrapping: usually 500-1000
    num_parallelCores = 6)                # total num of cores to use (parallelization)


# Extract new signatures (de-novo mutations)- may take a while
post.analysis <- 
  decipherMutationalProcesses(input = post.counts,
                              params = post.params)



#Downstream analyses and visualization
## examine the results
# Retrieve signatures (results)
post.sig <- post.analysis$Results$signatures

# Retrieve exposures (results)
post.exp <- post.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(post.sig, signature = 1, ylim = c(0, 0.10))
msigPlot(post.sig, signature = 2, ylim = c(0, 0.10))
msigPlot(post.sig, signature = 3, ylim = c(0, 0.10))
msigPlot(post.sig, signature = 4, ylim = c(0, 0.10))
#msigPlot(pre.sig, signature = 5, ylim = c(0, 0.10))




msigPlot(post.exp) + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))


# Export Signatures as data.frame
xprt <- coerceObj(x = post.sig, to = "data.frame") 
write.table(xprt, file = "")

# Get signatures from data (imported as data.frame) 
# and then convert it to mutSignatures object


# Retrieve COSMIC signatures from online repo, and then subset
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
cosmx<-cosmix[c(1, 2, 5,13)]


# match OVcar and COSMIC signatures
mSign.sar <- matchSignatures(mutSign = post.sig, reference = cosmix)
print(mSign.sar$plot)


#########
#pre.sig -> blcmx
#cosmix

blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = post.sig)




blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.2x <- blca.expo2$Results$count.result

# Plot exposures
bp1 <- msigPlot(blca.exp.1x, main = "BLCA | COSMIC sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

bp2 <- msigPlot(blca.exp.2x, main = "BLCA | pre-blca sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

# Visualize
grid.arrange(bp1, bp2, ncol = 2)




#####



qq<-data.frame(t(blca.exp.2x_df))

aa<-data.frame(qq[,1])
aa$sample_id<-rownames(qq)
aa$sig<-rep("sig.1")
colnames(aa)<-c("count","sample_id","sig")


bb<-data.frame(qq[,2])
bb$sample_id<-rownames(qq)
bb$sig<-rep("sig.2")
colnames(bb)<-c("count","sample_id","sig")


cc<-data.frame(qq[,3])
cc$sample_id<-rownames(qq)
cc$sig<-rep("sig.3")
colnames(cc)<-c("count","sample_id","sig")



dd<-data.frame(qq[,4])
dd$sample_id<-rownames(qq)
dd$sig<-rep("sig.4")
colnames(dd)<-c("count","sample_id","sig")

all<-rbind(aa,bb,cc,dd)

pdf("~/desktop/cexp_counts_annotated.pdf")
plota<-ggplot(all, aes(x = sample_id, y= count,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   



####




vv<-data.frame(t(xprt))

aa<-data.frame(vv[,1])
aa$sample_id<-rownames(vv)
aa$sig<-rep("sig.1")
colnames(aa)<-c("count","sample_id","sig")


bb<-data.frame(vv[,2])
bb$sample_id<-rownames(vv)
bb$sig<-rep("sig.2")
colnames(bb)<-c("count","sample_id","sig")


cc<-data.frame(vv[,3])
cc$sample_id<-rownames(vv)
cc$sig<-rep("sig.3")
colnames(cc)<-c("count","sample_id","sig")



dd<-data.frame(vv[,4])
dd$sample_id<-rownames(vv)
dd$sig<-rep("sig.4")
colnames(dd)<-c("count","sample_id","sig")

all_me<-rbind(aa,bb,cc,dd)

pdf("~/desktop/cexp_counts_annotated.pdf")
plotv<-ggplot(all_me, aes(x = sample_id, y= count,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plotv)
dev.off()   