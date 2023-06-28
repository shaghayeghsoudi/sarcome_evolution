#rm(list = ls())

## Title: Mutational signature analysis with mutSignature package in R (merged regions per patient)
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 
# Date: November 2022
# Updated: June 2023
## for vignette see: # https://cran.r-project.org/web/packages/mutSignatures/vignettes/get_sarted_with_mutSignatures.html

rm(list = ls())
library(rjson)
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)
library(reshape2)
library(kableExtra)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(ggpubr)
library(broom)
library(purrr)
source("http://peterhaschke.com/Code/multiplot.R")
# Load mutSignatures
library(mutSignatures)
# prep hg19
hg19 <- BSgenome.Hsapiens.UCSC.hg19

setwd("/Users/shsoudi/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/MutationalSignature")

#################################
### read and adjust meta table ##
#################################
meta_raw<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")

meta_raw$RTstatus<-ifelse(meta_raw$sequenceofsamplevRT== "beforeRT", "noRT",
                        ifelse(meta_raw$sequenceofsamplevRT== "afterRT", "postRT",
                        ifelse(meta_raw$sequenceofsamplevRT== "nopreopRT", "noRT",
                        "-")))

meta<-meta_raw%>%
  rename(sampleid="meta_id") %>%
  mutate(unique_sample_id=gsub("_.*$","",meta_id)) %>%
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

montecarlo_tumour_meta<-aa_montecarlo %>% 
  separate (SAMPLE,c("sample","sample_id","general_info")) %>%
  filter(grepl("T",general_info)) %>% ### subset muttaions only to tumour mutations
  mutate(general_info = gsub("T","",general_info), meta_id=paste(sample_id,general_info, sep = "_")) %>% 
  mutate(AF=TUMOR_DEPTH/TOTAL_DEPTH) %>%
  filter(AF > 0.05) %>%
  inner_join(meta,by = "meta_id") %>%
  mutate(chrom_pos=paste(CHR,POSITION , sep = "_"))%>%
  mutate(monto_ID=paste(identifier,chrom_pos , sep = "_")) %>%
  filter(!duplicated(monto_ID)) %>%
  select(CHR,POSITION,REF_ALLELE,TUMOR_ALLELE,identifier)%>%
  rename(CHR="CHROM" ,POSITION= "POS",REF_ALLELE="REF" ,TUMOR_ALLELE= "ALT",identifier="SAMPLEID")

## save output table
write.table(montecarlo_tumour_meta, file = "output_mutsignature_merged_regions_from_updated_solution_filterd_VAF0.05_montocarlo_June2023_tab.table", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


input_sigmutation<-read.table(file = "output_mutsignature_merged_regions_from_updated_solution_filterd_VAF0.05_montocarlo_June2023_tab.table", header = TRUE)
### filter out any muttaion that is not a SNV

muts<-c("A","T","C","G")
input_sigmutation_snv<-input_sigmutation[input_sigmutation$REF%in%muts,]
input_sigmutation_snv<-input_sigmutation_snv[input_sigmutation_snv$ALT%in%muts,]

##################################
#### start with Mutsignatures ####
##################################
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

write.table(y_pre, file = "outres_compute_mutType_MutSignatures_all_samples_merged_regions_VAF0.05.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

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

num.sign <- 2

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
pre.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 500,               # bootstrapping: usually 500-1000
    num_parallelCores = 10)                # total num of cores to use (parallelization)


# Extract new signatures (de-novo mutations)- may take a while
pre.analysis <- 
  decipherMutationalProcesses(input = pre.counts,
                              params = pre.params)


pre.sig <- pre.analysis$Results$signatures

# Retrieve exposures (results)
pre.exp <- pre.analysis$Results$exposures


# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
for (pp in 1:length(num.sign)){

  png(file = paste("denovo_MutSignature_sign_",pp,"_merged_regions.png", sep = ""))
  msigPlot(pre.sig, signature = pp, ylim = c(0, 0.10))
  dev.off()
}


# Export Signatures as data.frame
xprt <- coerceObj(x = pre.exp, to = "data.frame") 
write.table(xprt, file = "out_put_mutsignatures_xprt_denove_2_signature_count_All_RT.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
########################################################
##### costom visualization of denovo exposures #########
denovo_exp<-data.frame(t(xprt))

our_res_exp<-NULL
for(kk in 1:nrow(denovo_exp)){

  focal_sample<-denovo_exp[kk,]
  focal_sample<-(focal_sample/rowSums(focal_sample))*100
  count_table<-data.frame(t(focal_sample))
  count_table$sig<-rownames(count_table)
  count_table$sample_id<-rownames(denovo_exp[kk,])
  colnames(count_table)<-c("relative_contribution","sig","sample_id")
  rownames(count_table)<-NULL
  our_res_exp<-rbind(count_table,our_res_exp)
}

write.table(our_res_exp, file = "out_put_mutsignatures_xprt_relative_contribution_percentage_denove2_signature_All_RT.txt")

pdf("barplot_stacked_all_samples_2denovo_relative_contribution.pdf", width = 6, height = 7)
plota<-ggplot(our_res_exp, aes(x = sample_id, y= relative_contribution,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   

### boxplots 
our_res_exp$RTstatus<-gsub(".*_","",our_res_exp$sample_id)
our_res_exp$RTstatus<-gsub("naive","noRT",our_res_exp$RTstatus)
our_res_exp$RTstatus<-gsub("RTtreatment","postRT",our_res_exp$RTstatus)

denovos<-unique(our_res_exp$sig)

my.plot.denovo <- vector(mode = "list", length = 2)  ### adjust based on the number of Cosmic signatures taken

for(dd in 1:length(denovos)){
  focal<-our_res_exp[our_res_exp$sig==denovos[dd],]

  plot_sig1<-ggplot(focal, aes(x=RTstatus, y=relative_contribution, fill = RTstatus)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
                labs(title=denovos[dd]) 
                #F1<-plot_sig1 + stat_compare_means()
                F1 <- plot_sig1 + stat_compare_means(method = "t.test")
                my.plot.denovo[[dd]] <- F1
}

pdf("boxplot_ttest_denovo2_relative_contribution_percent_all_samples.pdf", width = 4, height = 5)
plots_denovo_all<-multiplot(plotlist = my.plot.denovo[1:2],cols= 1) 
print(plots_denovo_all)
dev.off()

#### heatmap of denovo exposure frequencies ####
freq<-data.frame(table(y_pre$SAMPLEID))

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
rownames(out_res_final_xprt)<-c("Sign.01","Sign.02")
out_res_final_xprt_df<-data.frame(out_res_final_xprt)

png("heatmap_denovo2_frequencies_all_samples.png", width = 700, height = 700)
pheatmap(out_res_final_xprt_df,cluster_rows =FALSE)
dev.off()


###############################################################
###############################################################
# Retrieve COSMIC signatures from online repo, and then subset
###############################################################
###############################################################
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
cosmx<-cosmix[c(1,2,5,13)]


# match OVcar and COSMIC signatures
png("heatmap_similarities_denovo_cosmic_all_samples.png", width = 700, height = 700)
mSign.sar <- matchSignatures(mutSign = pre.sig, reference = cosmix)
print(mSign.sar$plot)
dev.off()


blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = pre.sig)

blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.1x_df <- coerceObj(x = blca.exp.1x, to = "data.frame") 
write.table(blca.exp.1x_df, file = "output_mutsignatures_blca.exp.1x_cosmic_4signatures_count_All_samples_txt")

#blca.exp.2x <- blca.expo2$Results$count.result
#blca.exp.2x_df <- coerceObj(x = blca.exp.2x, to = "data.frame") 
#write.table(blca.exp.2x_df, file = "MutSignatures/Feb_2023/out_put_mutsignatures_blca.exp.2x_denove_signature_count_All_RT.txt")

# Plot exposures
cosmic_exp<-data.frame(t(blca.exp.1x_df))

our_res_exp<-NULL
for(kk in 1:nrow(cosmic_exp)){

  focal_sample<-cosmic_exp[kk,]
  focal_sample<-(focal_sample/rowSums(focal_sample))*100
  count_table<-data.frame(t(focal_sample))
  count_table$sig<-rownames(count_table)
  count_table$sample_id<-rownames(cosmic_exp[kk,])
  colnames(count_table)<-c("relative_contribution","sig","sample_id")
  rownames(count_table)<-NULL
  our_res_exp<-rbind(count_table,our_res_exp)
}

write.table(our_res_exp, file = "output_mutsignatures_normalized_blca.exp.1x_cosmic_4signatures_relative_contribution_All_samples_txt")


pdf("barplot_stacked_all_samples_4Cosmic_relative_contribution_all_samples.pdf", width = 6, height = 7)
plota<-ggplot(our_res_exp, aes(x = sample_id, y= relative_contribution,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   

### boxplots 
our_res_exp$RTstatus<-gsub(".*_","",our_res_exp$sample_id)
our_res_exp$RTstatus<-gsub("naive","noRT",our_res_exp$RTstatus)
our_res_exp$RTstatus<-gsub("RTtreatment","postRT",our_res_exp$RTstatus)
names(our_res_exp)[1]<-"percent_relative_contribution"

sigs<-unique(our_res_exp$sig)
my.plot <- vector(mode = "list", length = 4)  ### adjust based on the number of Cosmic signatures taken

for(ss in 1:length(sigs)){
  focal<-our_res_exp[our_res_exp$sig==sigs[ss],]

  plot_sig1<-ggplot(focal, aes(x=RTstatus, y=percent_relative_contribution, fill = RTstatus)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
                labs(title=sigs[ss]) 
                #F1<-plot_sig1 + stat_compare_means()
                F1 <- plot_sig1 + stat_compare_means(method = "t.test")
                my.plot[[ss]] <- F1
}

pdf("boxplot_ttest_Cosmic4_relative_contribution_percent_all_samples.pdf", width = 9, height = 10)
plots_cosmic_all<-multiplot(plotlist = my.plot[1:4],cols= 2) 
print(plots_cosmic_all)
dev.off()



#pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_cosmic_signatures_bp1_counts_all_samples.pdf")
#bp1 <- msigPlot(blca.exp.1x, main = "BLCA | COSMIC sigs.") + 
#  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))
#print(bp1)
#dev.off()


#pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_cosmic_signatures_bp2_counts_all_samples.pdf")
#bp2 <- msigPlot(blca.exp.2x, main = "BLCA | pre-blca sigs.") + 
#  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))
#print(bp1)
#dev.off()

# Visualize
#grid.arrange(bp1, bp2, ncol = 2)

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
#rownames(out_res_final_xprtcos)<-c("Cosmic.01","Cosmic.02","Cosmic.03","Cosmic.05","Cosmic.06","Cosmic.13")

out_res_final_xprtcos_df<-data.frame(out_res_final_xprtcos)

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/heatmap_4Cosmic_signatures_frequencies_all_samples_annotated_clustered.pdf", width = 8, height = 8)
plot_heat<-pheatmap(out_res_final_xprtcos_df,cluster_rows =FALSE)
print(plot_heat)
dev.off()

#ordered<-c("SRC125_preRT","SRC127_preRT","SRC167_preRT","SRC169_preRT","SRC170_preRT","SRC171_preRT","TB9051_preRT","SRC125_postRT","SRC127_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB22446_postRT","TB9051_postRT","RC130_noRT","SRC168_noRT" ,"SRC172_noRT","SRC173_noRT", "TB13092_noRT","TB13712_noRT","TB13959_noRT","TB9573_noRT")
#out_res_final_xprt_ordered<-out_res_final_xprt_df[order(match(colnames(out_res_final_xprt_df), ordered))]

#pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_ordered_Unclustered.pdf", width = 5, height = 5)
#pheatmap(out_res_final_xprt_ordered,cluster_cols=FALSE)
#dev.off()


##############
#### END #####
#############
