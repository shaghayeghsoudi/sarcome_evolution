

## Title: Mutational signature analysis with Muttaionalpatterns package in R
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 
# Date: November 2022
## for vignette see: https://bioconductor.org/packages/release/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.html

#rm(list = ls())
library(MutationalPatterns)
library(BSgenome)
library(ggplot2)
#library(dplyr)
library(stringr)
#library(tidyr)
library(ggpubr)
#head(available.genomes())

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")

meta<-read.delim(file = "metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE)
meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)


#### make required vcf file format
filenames_CN <- list.files("filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=TRUE, sep = "\t")[,c("Chromosome","Start","Ref", "Alt","AF","TumorReads","TumorDepth")]
     })

for (i in 1:length(attackStats_CN)){
    attackStats_CN[[i]]<-cbind(attackStats_CN[[i]],filenames_CN[i])
    }

### add sample name as a column
shared_clonal_mutations_fix <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix$sample_id<-gsub("filtered_variants/", "",gsub(".filtered.variants.oxomerge.final.txt","",shared_clonal_mutations_fix[,8]))
shared_clonal_mutations_fix$pos_id<-paste(shared_clonal_mutations_fix$Chromosome,shared_clonal_mutations_fix$Start, sep = "_")
shared_clonal_mutations_fix<-shared_clonal_mutations_fix[shared_clonal_mutations_fix$Start!="777428",]
shared_clonal_mutations_fix$unique_sample_id<-gsub("_.*$","",shared_clonal_mutations_fix$sample_id)
mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sampleid")


mutaion_RT$RTtretment<-ifelse(mutaion_RT$sequenceofsamplevRT== "beforeRT", "naive",
                        ifelse(mutaion_RT$sequenceofsamplevRT== "nopreopRT", "naive",
                        ifelse(mutaion_RT$sequenceofsamplevRT== "afterRT", "irradiated",
                        "-")))

mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id.x,mutaion_RT$RTtretment, sep = "_")
samples<-unique(mutaion_RT$identifier)

out_res<-NULL
for (i in 1:length(samples)){

    mutaion_RT_focal<-mutaion_RT[mutaion_RT$identifier==samples[i],]
    mutaion_RT_focal_nodup<-mutaion_RT_focal[!(duplicated(mutaion_RT_focal$pos_id)),]
    out_res<-rbind(mutaion_RT_focal_nodup,out_res)  
}

## convert to vcf file format
out_res$ID<-rep(".",nrow(out_res))
out_res$QUAL<-rep(".",nrow(out_res))
out_res$FILTER<-rep(".",nrow(out_res))
out_res$FORMAT<-rep("AD:DP",nrow(out_res))
out_res$Ref_count<-(out_res$TumorDepth-out_res$TumorReads)
out_res$allele<-paste(".",(paste(out_res$TumorReads,out_res$Ref_count, sep =":")),sep=",")
out_res$Ref_count<-paste("t_ref_count=",out_res$Ref_count,sep = "")
out_res$TumorReads<-paste("t_alt_count=",out_res$TumorReads,sep = "")
out_res$INFO<-paste(out_res$TumorReads,out_res$Ref_count, sep =";" )
#out_res$Chromosome<-gsub("chr","",out_res$Chromosome)
#muts<-c("A","T","C","G")
#out_res<-out_res[out_res$Ref%in%muts,]
#out_res<-out_res[out_res$Alt%in%muts,]

#out_res$Chromosome<-gsub("chr", "",out_res$Chromosome)
out_res$id1<-paste(out_res$Chromosome,out_res$Start, sep = ":")
out_res$id2<-paste(out_res$Ref,out_res$Alt, sep = "/")
out_res$ID<-paste(out_res$id1,out_res$id2, sep = "_")


muttaion_selected_cols<-out_res[,c("Chromosome","Start","ID","Ref","Alt","QUAL","FILTER","INFO","FORMAT","allele","identifier")]
colnames(muttaion_selected_cols)<-c("#CHROM","POS","ID","REF","ALT"	,"QUAL"	,"FILTER",	"INFO",	"FORMAT","SAMPLE", "sample_id")
colnames(muttaion_selected_cols)<-c("#CHROM","POS","ID","REF","ALT"	,"QUAL"	,"FILTER",	"INFO",	"FORMAT","SAMPLE", "sample_id")

samples<-unique(muttaion_selected_cols$sample_id)


for (ii in 1:length(samples)){

    focal_vcf<-muttaion_selected_cols[muttaion_selected_cols$sample_id==samples[ii],]
    colnames(focal_vcf)[10]<-samples[ii]
    focal_vcf<-focal_vcf[,-11]
    #focal_vcf<- focal_vcf[!grepl("-",focal_vcf$REF),]
    #focal_vcf<- focal_vcf[!grepl("-",focal_vcf$ALT),]
    write.table(focal_vcf, file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/vcf_format_with_indels/",samples[ii],".filtered.variants.oxomerge.final.Jan2023.vcf", sep = ""),row.names= FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

##fileformat=VCFv4.1
#for (v in length(samples)){
#prettyprint <- function() {
#  cat('##fileformat=VCFv4.1',
#  '\n##FORMAT=<ID=AD,Number=2,Type=Integer,Description= "Allelic depths (number of reads in each observed allele)" >',
#  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth" >',
#  '\n##FORMAT=<ID=FT,Number=1,Type=String,Description="Variant filters" >',
#  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Tumour ref count" >',
#  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Tumour alt counth" > \n'  
#      )
#      print(muttaion_selected_cols$sample_id==samples[v])
# }
# prettyprint()
#}


### example how vcf format should look like
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PD21928b
#1	105605108	1:105605108_T/A	T	A	429.75	PASS	.	GT:AD:DP:GQ:PL	0/0:29,0:29:75:0,75,818
#1	120042428	1:120042428_T/C	T	C	440.79	PASS	.	GT:AD:DP:GQ:PL	0/0:27,0:27:78:0,78,862


####################################
#### run the muttaionalpatterns ####
vcf_files <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/vcf_format_with_indels", pattern = "*.filtered.variants.oxomerge.final.Jan2023.vcf",full.names = TRUE)
sample_names <- samples
sar_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, predefined_dbs_mbs = TRUE)  ## Any neighbouring SNVs will be merged into DBS/MBS variants.
                                                                      ## Set the 'predefined_dbs_mbs' to 'TRUE' if you don't want this.
#sar_snv_grl <- get_mut_type(sar_grl, type = "snv")
#sar_indel_grl <- get_mut_type(sar_grl, type = "indel")
#sar_dbs_grl <- get_mut_type(sar_grl, type = "dbs")
#sar_mbs_grl <- get_mut_type(sar_grl, type = "mbs")

#indel_grl <- read_vcfs_as_granges(vcf_files, sample_names, 
#                              ref_genome, type = "indel")

######
##indel_grl <- get_indel_context(indel_grl, ref_genome)
#head(indel_grl[[1]], n = 3)

#Next count the number of indels per type. This results in a matrix that is similar to the mut_mat matrix.

#indel_counts <- count_indel_contexts(indel_grl)
#head(indel_counts)
#plot the indel spectra
#plot_indel_contexts(indel_counts, condensed = TRUE)

#pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/Mutation_spectrum_exome_indels_Feb2023.pdf", width = 10, height = 30)
#indels<-plot_main_indel_contexts(indel_counts)
#print(indels)
#dev.off()
#Using other signature matrixes

#####
muts <- mutations_from_vcf(sar_grl[[1]])
head(muts, 12)

types <- mut_type(sar_grl[[1]])
head(types, 12)

context <- mut_context(sar_grl[[1]], ref_genome)
head(context, 12)

type_context <- type_context(sar_grl[[1]], ref_genome)
lapply(type_context, head, 12)

type_occurrences <- mut_type_occurrences(sar_grl, ref_genome)
write.table(type_occurrences, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/type_occurrences_muutaionpattern_exome_noindels_Feb2023.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


## customized plots
type_occurrences$RTstatus<-gsub(".*_","",rownames(type_occurrences))
rownames(type_occurrences)<-NULL

types<-c("C>A","C>G","C>T","T>A","T>C","T>G","C>T at CpG","C>T other")
out_res<-NULL
for (ii in 1:length(types)){
    focal_type<-data.frame("count"=type_occurrences[,ii])
    focal_type$type<-rep(names(type_occurrences[ii]))
    focal_type$RTstatus<-type_occurrences$RTstatus
    out_res<-rbind(focal_type,out_res) 

  }

out_res$RTstatus<-gsub("preRT","noRT",out_res$RTstatus)
box_type<-ggplot(out_res, aes(x=type, y=count, fill=RTstatus)) + 
    geom_boxplot()
ggsave(filename="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/summary_all_WES_noindels_regions_added_Feb2023/boxplot_nucleotide_types_per_radiation.jpg", plot=box_type, height = 5, width = 8)


## 96 muttaional profile
mut_mat <- mut_matrix(vcf_list = sar_grl, ref_genome = ref_genome)
head(mut_mat)
plot_96_profile(mut_mat[, c(1: 7)])


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/Mutation_spectrum_exome_96profile_noindel_Jan2023.pdf", width = 7, height = 15)
plot_A<-plot_96_profile(mut_mat[, c(1:27)])
print(plot_A)
dev.off()

#### De novo mutational signature extraction using NMF
#### denovo signature
# 1: NMF
mut_mat <- mut_mat + 0.0001
library("NMF")
estimate <- nmf(mut_mat, rank = 2:6, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/NMF_optimal_factorization_rank_exome_96profile_noindel_Feb2023.pdf", width = 5, height = 5)
estimate<-plot(estimate)
plot(estimate)
print(estimate)
dev.off()

### I chose three signature 
### the most common approach is to choose the smallest rank for which cophenetic correlation coefficient starts decreasing
nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 20, single_core = TRUE)
#rownames(nmf_res$contribution)<-c("SBS1-like","SBSA","SBSB")
#colnames(nmf_res$signatures) <- c("SBS1-like","SBSA","SBSB")


#It’s possible that some of the signatures extracted by NMF are very similar to signatures 
#that are already known. Therefore, it might be useful to change the names of the NMF 
#signatures to these already known signatures. This often makes it easier to interpret 
#your results.
signatures = get_known_signatures()
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
#[1] "SBSA"      "SBS1-like" "SBSB" 

nmf_res_contribution<-data.frame(nmf_res$contribution)
#rownames(nmf_res_contribution)<-c("SBS1-like","SBSA","SBSB")
write.table(nmf_res_contribution, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/out_res_nmf_res_contribution_exome_96profile_noindel_Feb2023.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)


out_res_nmf<-NULL
for (jj in 1:ncol(nmf_res_contribution)){

  nmf_focal<-data.frame("contribution"= nmf_res_contribution[,jj])
  normalized<-1/colSums(nmf_focal)
  nmf_focal$relative_contribution<-(nmf_focal$contribution)*normalized
  nmf_focal$signature<-rownames(nmf_res_contribution)
  nmf_focal$patient_ID<-colnames(nmf_res_contribution[jj])
  out_res_nmf<-rbind( nmf_focal,out_res_nmf)

}
write.table(out_res_nmf, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/out_res_normalized_nmf_res_contribution_exome_96profile_noindel_Feb2023.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)


png(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/96profile_noindel_Feb2023.png", width = 700, height = 700)
plot_96profile<-plot_96_profile(nmf_res$signatures, condensed = TRUE)
print(plot_96profile)
dev.off()


## make bar plot of relative contribution of each signature
out_res_nmf$RTstatus<-gsub(".*_","",out_res_nmf$patient_ID)
denovo_contribution<-ggplot(out_res_nmf, aes(fill=signature , y=relative_contribution, x=patient_ID)) + 
    geom_bar(position="fill", stat="identity")
     denovo_contribution + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(filename="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/barplot_mutpatterns_relative_contribution_NMF_denovo_SNVs.png", height = 500, width = 500)


##### compare if there is any difference between signatures among different pre,post and noRT samples ###
#out_res_nmf$RTstatus<-gsub("preRT","noRT",out_res_nmf$RTstatus)
SBS1.like<-out_res_nmf[out_res_nmf$signature=="SBS1-like",]
SBSA<-out_res_nmf[out_res_nmf$signature=="SBSA",]
SBSB<-out_res_nmf[out_res_nmf$signature=="SBSB",]



### box plots of difference betwen groups in each signature ###
#res <- wilcox.test(relative_contribution~ RTstatus,
#                   data = SBS1.like,
#                   exact = FALSE)

plot_sig1<-ggplot(SBS1.like, aes(x=RTstatus, y=relative_contribution, fill = RTstatus)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
                labs(title="SBS1.like") 
                F1<-plot_sig1 + stat_compare_means()
                #F1 <- plot_sig1 + stat_compare_means(method = "t.test")



plot_sig2<-ggplot(SBSA, aes(x=RTstatus, y=relative_contribution, fill = RTstatus)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)+
                labs(title="SBSA")
                F2<-plot_sig2 + stat_compare_means()


plot_sig3<-ggplot(SBSB, aes(x=RTstatus, y=relative_contribution, fill = RTstatus)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)+
                labs(title="SBSB")
                F3<-plot_sig3 + stat_compare_means()


### check normality
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/histogram_denovo_relative_contribution.pdf", width = 7, height = 5)
par(mfrow=c(1,3)) 
hist(SBS1.like$relative_contribution, breaks = 10, main = "SBS1.lik")
hist(SBSA$relative_contribution, breaks = 10, main = "SBSA")
hist(SBSB$relative_contribution, breaks = 10, main = "SBSB")
dev.off()
shapiro.test(SBSB$relative_contribution)




library("gridExtra")
png("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/boxplot_relative_contribution_customized_WED_SNVs_Feb2023.png", width = 700, height = 1000, res = 150)
plot_all_box<-grid.arrange(F1,F2,F3, nrow= 3)
print(plot_all_box)
dev.off()


#combi_mat = rbind(indel_counts, dbs_counts)
#nmf_res_combi <- extract_signatures(combi_mat, rank = 2, nrun = 10, single_core = TRUE)

#B: Bayesian NMF
# BiocManager::install("ccfindR")

### get similar signatures with Cosmic
signatures = get_known_signatures()
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)

#Visualizing NMF results

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/96profile_noindel_NMF_exome_Jan2023.pdf", width = 7, height = 7)
sig_profile<-plot_96_profile(nmf_res$signatures, condensed = TRUE)
print(sig_profile)
dev.off()


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/barplot_contribution_noindel_NMF_exome_Jan2023.pdf", width = 5, height = 5)
contribution<-plot_contribution(nmf_res$contribution, nmf_res$signature,
  mode = "relative"
)
print(contribution)
dev.off()


### plot the heatmap

hclust_signatures <- cluster_signatures(nmf_res$signatures, method = "average")
signatures_order <- colnames(nmf_res$signatures)[hclust_signatures$order]
signatures_order


hclust_samples <- cluster_signatures(mut_mat, method = "average")
samples_order <- colnames(mut_mat)[hclust_samples$order]
samples_order


pdf(file = "~Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/heatmap_contribution_noindel_NMF_exome_unclustered.pdf", width = 7, height = 5)
heatmap_contribution<-plot_contribution_heatmap(nmf_res$contribution,
  sig_order = signatures_order, sample_order = samples_order,
  cluster_sigs = FALSE, cluster_samples = FALSE
)
print(heatmap_contribution)
dev.off()


plot_compare_profiles(mut_mat[, 1],
  nmf_res$reconstructed[, 1],
  profile_names = c("Original", "Reconstructed"),
  condensed = TRUE
)


#########################
## Signature refitting
#Find mathematically optimal contribution of COSMIC signatures
fit_res <- fit_to_signatures(mut_mat, signatures)
fit_contribution<-data.frame(fit_res$contribution)
write.table(fit_contribution, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/out_res_fit_to_cosmic_contribution_signature_WES_noindel_Jan2023.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/barplot_fit_to_cosmic_contribution_signature_WES_noindel_Jan2023.pdf", width = 10, height = 10)
cos_contribution<-plot_contribution(fit_res$contribution,
  coord_flip = FALSE,
  mode = "absolute"
)
print(cos_contribution)
dev.off()


samples<-colnames(fit_contribution)
#samples<-samples[c(-6,-18)]
fit_contribution$cosmic<-rownames(fit_contribution)

out_res_cos<-NULL
for(ss in 1:length(samples)){

   fit_contribution_focal<-fit_contribution[,c(samples[ss],"cosmic")]
   colnames(fit_contribution_focal)<-c("contribution","cosmic")
   fit_contribution_focal_order<-fit_contribution_focal[order(fit_contribution_focal$contribution,decreasing = TRUE),]   
   fit_contribution_focal_order<-fit_contribution_focal_order[fit_contribution_focal_order$contribution>0,]
   fit_contribution_focal_order<-fit_contribution_focal_order[1:5,]
   fit_contribution_focal_order$sample_id<-rep(samples[ss])
   out_res_cos<-rbind(out_res_cos,fit_contribution_focal_order)
}

rownames(out_res_cos)<-NULL
out_res_cos<-out_res_cos[out_res_cos$sample_id!="cosmic",]
out_res_cos$contribution<-as.integer(out_res_cos$contribution)

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/barplot_cosmic_contribution_customized_Jan2023.pdf", width = 5, height = 5)
plotcosmic<-ggplot(out_res_cos, aes(x = sample_id, y= contribution,fill = cosmic)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plotcosmic)
dev.off() 


#### box plot of cosmic differences bteween pre and postRT samples 
out_res_cos2<-NULL
for(ss in 1:length(samples)){

   fit_contribution_focal<-fit_contribution[,c(samples[ss],"cosmic")]
   colnames(fit_contribution_focal)<-c("contribution","cosmic")
   fit_contribution_focal_order<-fit_contribution_focal[order(fit_contribution_focal$contribution,decreasing = TRUE),]   
   fit_contribution_focal_order<-fit_contribution_focal_order[fit_contribution_focal_order$contribution>=0,]
   #fit_contribution_focal_order<-fit_contribution_focal_order[1:5,]
   fit_contribution_focal_order$sample_id<-rep(samples[ss])
   out_res_cos2<-rbind(out_res_cos2,fit_contribution_focal_order)
}

rownames(out_res_cos2)<-NULL
out_res_cos2<-out_res_cos2[out_res_cos2$sample_id!="cosmic",]
out_res_cos2$contribution<-as.integer(out_res_cos2$contribution)

out_res_cos2_pre_post<-out_res_cos2[grepl("pre",out_res_cos2$sample_id)|grepl("post",out_res_cos2$sample_id),]
out_res_cos2_pre_post$RT_status<-gsub(".*_","",out_res_cos2_pre_post$sample_id)

sbs1<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS1",]
sbs13<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS13",]
sbs2<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS2",]
sbs3<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS3",]
sbs5<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS5",]
sbs8<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS8",]
sbs7a<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS7a",]
sbs15<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS15",]

cosmic1<-ggplot(sbs1, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="refitted cosmic SBS1")
                cosmic1_stat<-cosmic1 + stat_compare_means()


cosmic2<-ggplot(sbs2, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="refitted cosmic SBS2")
                cosmic2_stat<-cosmic2 + stat_compare_means()

cosmic3<-ggplot(sbs3, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="refitted cosmic SBS3")
                cosmic3_stat<-cosmic3 + stat_compare_means()                


cosmic5<-ggplot(sbs5, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="refitted cosmic SBS5")
                cosmic5_stat<-cosmic5 + stat_compare_means()       


cosmic8<-ggplot(sbs8, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="refitted cosmic SBS8")
                cosmic8_stat<-cosmic8 + stat_compare_means()     


cosmic13<-ggplot(sbs13, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="refitted cosmic SBS13")
                cosmic13_stat<-cosmic13 + stat_compare_means()   


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/boxplot_cosmic_contribution_between_pre_post_compsrison_Jan2023.pdf", width = 8, height = 12)
plot_cosmics<-grid.arrange(cosmic1_stat,cosmic2_stat,cosmic3_stat,cosmic5_stat,cosmic8_stat,cosmic13_stat,ncol=2)
print(plot_cosmics)
dev.off()




##############################
## Stricter refitting
### Bootstrapped refitting.

contri_boots <- fit_to_signatures_bootstrapped(mut_mat[, c(3, 7)],
  signatures,
  n_boots = 50,
  method = "strict"
)
plot_bootstrapped_contribution(contri_boots)


plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "dotplot")



 #Similarity between mutational profiles and signatures         

cos_sim(mut_mat[, 1], signatures[, 1])
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, signatures)
cos_sim_samples_signatures[1:3, 1:3]                    