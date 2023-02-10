#rm(list = ls())

## Title: Mutational signature analysis with mutSignature package in R (merged regions per patient)
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 
# Date: November 2022
# Updated: February 2023
## for vignette see: # https://cran.r-project.org/web/packages/mutSignatures/vignettes/get_sarted_with_mutSignatures.html


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

### make input file for mutsignatures ###
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
meta<-read.delim(file = "metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE)
meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)

mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sampleid")
mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]


mutaion_RT$RTtretment<-ifelse(mutaion_RT$sequenceofsamplevRT== "beforeRT", "naive",
                        ifelse(mutaion_RT$sequenceofsamplevRT== "nopreopRT", "naive",
                        ifelse(mutaion_RT$sequenceofsamplevRT== "afterRT", "RTtreatment",
                        "-")))

mutaion_RT<-mutaion_RT[mutaion_RT$Chromosome!="Start",]
#cols <- c("Impact","Coding")
#mutaion_RT[cols] <- lapply(mutaion_RT[cols], as.character)                   
#mutaion_RT$Impact[mutaion_RT$Coding== "splicing"] <- "Splice sites"
#colnames(mutaion_RT)<-c("sampleID","chr","pos","ref","mut","RT_code","Gene","Impact","RTtretment","RT_staus_code")

mutaion_RT$pos_id<-paste(mutaion_RT$Chromosome,mutaion_RT$Start, sep = "_")
mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id,mutaion_RT$RTtretment, sep = "_")


samples<-unique(mutaion_RT$identifier)

out_res<-NULL
for (i in 1:length(samples)){

    mutaion_RT_focal<-mutaion_RT[mutaion_RT$identifier==samples[i],]
    mutaion_RT_focal_nodup<-mutaion_RT_focal[!(duplicated(mutaion_RT_focal$pos_id)),]
    out_res<-rbind(mutaion_RT_focal_nodup,out_res)  
}

out_res_good<-out_res[,c("Chromosome", "Start","Ref","Alt","identifier")]
colnames(out_res_good)<-c("CHROM","POS","REF","ALT","SAMPLEID")
out_res_good$CHROM<-paste("chr",out_res_good$CHROM,sep = "")

## write input file for mutsignatures ###
write.table(out_res_good, file = "MutSignatures/Feb_2023/outres_input_mutsignatures_merged_regions_all_samples_from_updated_solution_metadata_Feb2023_tab.table", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



input_sigmutation<-read.table(file = "MutSignatures/Feb_2023/outres_input_mutsignatures_merged_regions_all_samples_from_updated_solution_metadata_Feb2023_tab.table", header = TRUE)
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
#y_pre<-input_sigmutation_snv[grepl("noRT",input_sigmutation_snv$SAMPLEID),]
y_pre<-input_sigmutation_snv
#De novo extraction of Mutational Signatures 
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



write.table(y_pre, file = "MutSignatures/Feb_2023/outres_compute_mutType_MutSignatures_all_samples_merged_regions.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

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

num.sign <- 3

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
## examine the results
# Retrieve signatures (results)
pre.sig <- pre.analysis$Results$signatures

# Retrieve exposures (results)
pre.exp <- pre.analysis$Results$exposures


# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
for (pp in 1:length(num.sign)){

  png(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/denovo_MutSignature_sign_",pp,"_PreRT_merged_regions_boot500_cgDNA.png"))
  msigPlot(pre.sig, signature = pp, ylim = c(0, 0.10))
  dev.off()
}


png(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_denovo_MutSignature_counts_per_sample_all_samples_merged_regions_boot500.png")
msigPlot(pre.exp) + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))
dev.off()

write.table(xprt, file = "MutSignatures/Feb_2023/out_put_mutsignatures_xprt_denove_3_signature_count_All_RT.txt")


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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_cexp_counts_all_samples_sample_annotated.png", width = 500, height = 500)
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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/heatmap_denovo_frequencies_all_samples_annotated_clustered.png", width = 700, height = 700)
pheatmap(out_res_final_xprt_df,cluster_rows =FALSE)
dev.off()

### make ordered heatmap ###
#ordered<-c("SRC125_preRT","SRC127_preRT","SRC167_preRT","SRC169_preRT","SRC170_preRT","SRC171_preRT","TB9051_preRT","SRC125_postRT","SRC127_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB22446_postRT","TB9051_postRT","RC130_noRT","SRC168_noRT" ,"SRC172_noRT","SRC173_noRT", "TB13092_noRT","TB13712_noRT","TB13959_noRT","TB9573_noRT")
#out_res_final_xprt_ordered<-out_res_final_xprt_df[order(match(colnames(out_res_final_xprt_df), ordered))]

#pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_ordered_Unclustered.pdf", width = 5, height = 5)
#pheatmap(out_res_final_xprt_ordered,cluster_cols=FALSE)
#dev.off()

for (jj in 1:nrow(out_res_final_xprt)){
  
  sig1_deno<-data.frame("proportion"=out_res_final_xprt[jj,])
  sig1_deno$sample_id<-rownames(sig1_deno)
  sig1_deno<-sig1_deno%>%separate(sample_id,c("uniq_sample_id","RT"))
  
  png(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/boxplot_denovo_MutSignature_sign_",jj,"_all_samples_merged_boot500_cgDNA.png"), width = 500, height = 500)
  plot_sig1<-ggplot(sig1_deno, aes(x=RT, y=proportion, fill = RT)) + 
  geom_boxplot(outlier.colour="black",
                 outlier.size=0.7) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)+
  labs(title=rownames(out_res_final_xprt)[jj])
  sig1<-plot_sig1 + stat_compare_means()
  print(sig1)
  dev.off()
}
     
###############################################################
###############################################################
# Retrieve COSMIC signatures from online repo, and then subset
###############################################################
###############################################################
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
cosmx<-cosmix[c(1, 2, 5,13)]

# match OVcar and COSMIC signatures
png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/heatmap_similarities_denovo_cosmic_all_samples.png", width = 700, height = 700)
mSign.sar <- matchSignatures(mutSign = pre.sig, reference = cosmix)
print(mSign.sar$plot)
dev.off()


blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = pre.sig)

blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.1x_df <- coerceObj(x = blca.exp.1x, to = "data.frame") 
write.table(blca.exp.1x_df, file = "MutSignatures/Feb_2023/out_put_mutsignatures_blca.exp.1x_cosmic_signature_count_All_RT.txt")


blca.exp.2x <- blca.expo2$Results$count.result
blca.exp.2x_df <- coerceObj(x = blca.exp.2x, to = "data.frame") 
write.table(blca.exp.2x_df, file = "MutSignatures/Feb_2023/out_put_mutsignatures_blca.exp.2x_denove_signature_count_All_RT.txt")


# Plot exposures
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_cosmic_signatures_bp1_counts_all_samples.pdf")
bp1 <- msigPlot(blca.exp.1x, main = "BLCA | COSMIC sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))
print(bp1)
dev.off()


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_cosmic_signatures_bp2_counts_all_samples.pdf")
bp2 <- msigPlot(blca.exp.2x, main = "BLCA | pre-blca sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))
print(bp1)
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

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_cosmic_signatures_counts_all_samples_annotated.png", width = 700, height = 700)
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
rownames(out_res_final_xprtcos)<-c("Cosmic.01","Cosmic.02","Cosmic.05","Cosmic13")
out_res_final_xprtcos_df<-data.frame(out_res_final_xprtcos)

png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/heatmap_Cosmic_signatures_frequencies_all_samples_annotated_clustered.png", width = 700, height = 700)
plot_heat<-pheatmap(out_res_final_xprtcos_df,cluster_rows =FALSE)
print(plot_heat)
dev.off()

#ordered<-c("SRC125_preRT","SRC127_preRT","SRC167_preRT","SRC169_preRT","SRC170_preRT","SRC171_preRT","TB9051_preRT","SRC125_postRT","SRC127_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB22446_postRT","TB9051_postRT","RC130_noRT","SRC168_noRT" ,"SRC172_noRT","SRC173_noRT", "TB13092_noRT","TB13712_noRT","TB13959_noRT","TB9573_noRT")
#out_res_final_xprt_ordered<-out_res_final_xprt_df[order(match(colnames(out_res_final_xprt_df), ordered))]

#pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_ordered_Unclustered.pdf", width = 5, height = 5)
#pheatmap(out_res_final_xprt_ordered,cluster_cols=FALSE)
#dev.off()

