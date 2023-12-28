


#rm(list = ls())
library(MutationalPatterns)
library(BSgenome)
library(dplyr)
library(stringr)
library(tidyr)
library("ggplot2")
library("gridExtra")
#head(available.genomes())

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)


setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated.txt", header = FALSE, sep= "\t")
colnames(meta)<-c("sample_id","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)


filenames_CN <- list.files("cfDNA_analysis/Updated_files_pyclone/pyclone_multiple_Titan_probesadd500_puritycutoff_fullsegs_final_full/tsv",pattern="*pyclone.results.tsv", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })


py_data <- do.call("rbind", attackStats_CN) 
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))
py_data_good_RTID<-merge(py_data_good,meta, by.x = "sample_id", by.y= "sample_id")
py_data_good_RTID$Chrom<-gsub("chr","",py_data_good_RTID$Chrom)
py_data_good_RTID$pos_id<-paste(py_data_good_RTID$Chrom,py_data_good_RTID$Pos,sep = ":")
py_data_good_RTID$allele_id<-paste(py_data_good_RTID$Ref,py_data_good_RTID$Alt, sep = "/")
py_data_good_RTID$ID<-rep(".",nrow(py_data_good_RTID))
py_data_good_RTID$QUAL<-rep(".",nrow(py_data_good_RTID))
py_data_good_RTID$FILTER<-rep(".",nrow(py_data_good_RTID))
py_data_good_RTID$FORMAT<-rep("AD:DP",nrow(py_data_good_RTID))
py_data_good_RTID$INFO<-rep("t_alt_count")
py_data_good_RTID$sample<-rep("..:..")

samples<-unique(py_data_good_RTID$unique_sample_id)

for(i in 1:length(samples)){

    py_focal<-py_data_good_RTID[py_data_good_RTID$unique_sample_id==samples[i],]
    py_focal$ID<-paste(py_focal$pos_id,py_focal$allele_id, sep = "_")
    #py_focal_ordered<-py_focal[order(py_focal$cellular_prevalence, decreasing = TRUE),]
    find_clonal<-py_focal[which.max(py_focal$cellular_prevalence),]
    trunck_muttaions<-py_focal[py_focal$cluster_id==find_clonal$cluster_id,]
    trunck_muttaions_uniq<-trunck_muttaions[!duplicated(trunck_muttaions$pos_id),]
    trunck_muttaions_uniq<-trunck_muttaions_uniq[,c("Chrom","Pos","ID", "Ref","Alt","QUAL","FILTER","INFO","FORMAT","sample")]
    colnames(trunck_muttaions_uniq)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample")
    colnames(trunck_muttaions_uniq)[10]<-samples[i]
    write.table(trunck_muttaions_uniq,file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/cfDNA_analysis/Updated_files_pyclone/vcfs_trunck_branch/vcfs_trunck_muttaions_",samples[i],"_cfDNA_pyclone.vcf", sep = ""),col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

    branch_muttaions<-py_focal[py_focal$cluster_id!=find_clonal$cluster_id,]
    branch_muttaions_uniq<-branch_muttaions[!duplicated(branch_muttaions$pos_id),]
    branch_muttaions_uniq<-branch_muttaions_uniq[,c("Chrom","Pos","ID", "Ref","Alt","QUAL","FILTER","INFO","FORMAT","sample")]
    colnames(branch_muttaions_uniq)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample")
    colnames(branch_muttaions_uniq)[10]<-samples[i]
    write.table(branch_muttaions_uniq,file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/cfDNA_analysis/Updated_files_pyclone/vcfs_trunck_branch/vcfs_branch_muttaions_",samples[i],"_cfDNA_pyclone.vcf", sep = ""),col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)



}


################################
### run mutationalPatterns ####


vcf_files <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/cfDNA_analysis/Updated_files_pyclone/vcfs_trunck_branch",pattern = "*_cfDNA_pyclone.vcf",full.names = TRUE)
sample_names <- c("SRC125_branch","SRC127_branch","SRC130_branch","SRC150_branch","SRC167_branch","SRC168_branch","SRC169_branch","SRC170_branch","SRC171_branch","SRC172_branch","SRC173_branch","TB12052_branch","TB13092_branch","TB13712_branch","TB13959_branch","TB22446_branch","TB8016_branch","TB9051_branch","TB9573_branch","SRC125_trunck","SRC127_trunck","SRC130_trunck","SRC150_trunck","SRC167_trunck","SRC168_trunck","SRC169_trunck","SRC170_trunck","SRC171_trunck","SRC172_trunck","SRC173_trunck","TB12052_trunck","TB13092_trunck","TB13712_trunck","TB13959_trunck","TB22446_trunck","TB8016_trunck","TB9051_trunck","TB9573_trunck")
sar_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, predefined_dbs_mbs = TRUE)


muts <- mutations_from_vcf(sar_grl[[1]])
types <- mut_type(sar_grl[[1]])
context <- mut_context(sar_grl[[1]], ref_genome)
type_context <- type_context(sar_grl[[1]], ref_genome)
lapply(type_context, head, 12)

type_occurrences <- mut_type_occurrences(sar_grl, ref_genome)
write.table(type_occurrences, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/pyclone_trunck_vs_branch_type_occurrences_muutaionpattern_cfDNA.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


type_occurrences_good<-type_occurrences[,1:6]
type_occurrences_good$pyclone_status<-gsub(".*_","",rownames(type_occurrences_good))
mut_types<-c("C>A","C>G","C>T","T>A","T>C","T>G")

out_mut<-NULL
for(mm in 1:length(mut_types)){

    mut_focal<-type_occurrences_good[,c(mm,7)]
    mut_focal$mut$type<-rep(mut_types[mm])
    colnames(mut_focal)<-c("mut_count","RT_status","mut_type")
    rownames(mut_focal)<-NULL 
    #mut_focal$mut_count<-as.numeric(mut_focal[,1])
    out_mut<-rbind(mut_focal,out_mut)
}

out_mut$mut_type<- as.character(out_mut$mut_type)
out_mut$mut_count<- as.integer(out_mut$mut_count)

plotA<-ggplot(data = out_mut, aes(x = mut_type, y = mut_count, fill = RT_status, color = RT_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))


out_mut_pre_post<-out_mut[grepl("pre",out_mut$RT_status) | grepl ("post",out_mut$RT_status),]

plotB<-ggplot(data = out_mut_pre_post, aes(x = mut_type, y = mut_count, fill = RT_status, color = RT_status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))


test_sig<-out_mut_pre_post[out_mut_pre_post$mut_type==mut_types[mm],]
res <- wilcox.test(mut_count~ RT_status,
                   data = test_sig,
                   exact = FALSE)


res
###


type_occurrences_good_p<-type_occurrences[,-3]
type_occurrences_good_p$RT_status<-gsub(".*_","",rownames(type_occurrences_good_p))
mut_types<-c("C>A", "C>G", "T>A", "T>C","T>G","C>T_at_CpG","C>T_other")

out_mut_p<-NULL
for(mm in 1:length(mut_types)){

    mut_focal<-type_occurrences_good_p[,c(mm,8)]
    mut_focal$mut$type<-rep(mut_types[mm])
    colnames(mut_focal)<-c("mut_count","RT_status","mut_type")
    rownames(mut_focal)<-NULL 
    #mut_focal$mut_count<-as.numeric(mut_focal[,1])
    out_mut_p<-rbind(mut_focal,out_mut_p)
}

out_mut_p$mut_type<- as.character(out_mut_p$mut_type)
out_mut_p$mut_count<- as.integer(out_mut_p$mut_count)

plotC<-ggplot(data = out_mut_p, aes(x = mut_type, y = mut_count, fill = RT_status, color = RT_status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))


out_mut_pre_post_p<-out_mut_p[grepl("pre",out_mut_p$RT_status) | grepl ("post",out_mut_p$RT_status),]

plotD<-ggplot(data = out_mut_pre_post_p, aes(x = mut_type, y = mut_count, fill = RT_status, color = RT_status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))


test_sig_p<-out_mut_pre_post_p[out_mut_pre_post_p$mut_type==mut_types[mm],]
res <- wilcox.test(mut_count~ RT_status,
                   data = test_sig_p,
                   exact = FALSE)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/boxplot_count_mutation_type_cfDNA.pdf", width = 10, height = 7)
plot_all<-grid.arrange(plotA,plotB,plotC,plotD, ncol = 2)
print(plot_all)
dev.off()


##########
p1 <- plot_spectrum(type_occurrences)
p2 <- plot_spectrum(type_occurrences, CT = TRUE)
p3 <- plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE, legend = FALSE)


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/Mutation_spectrum_exome_all_RT_cfDNA.pdf", width = 12, height = 5)
plot_all<-grid.arrange(p1, p2, p3, ncol = 3)
print(plot_all)
dev.off()


pre<-type_occurrences[grepl("pre", rownames(type_occurrences)),]
p1_pre <- plot_spectrum(pre)
p2_pre <- plot_spectrum(pre, CT = TRUE)
p3_pre <- plot_spectrum(pre, CT = TRUE, 
                    indv_points = TRUE, legend = FALSE)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/Mutation_spectrum_cfDNA_preRT.pdf", width = 12, height = 5)
plot_pre<-grid.arrange(p1_pre, p2_pre, p3_pre, ncol = 3)
print(plot_pre)
dev.off()


post<-type_occurrences[grepl("post", rownames(type_occurrences)),]
p1_post <- plot_spectrum(post)
p2_post <- plot_spectrum(post, CT = TRUE)
p3_post <- plot_spectrum(post, CT = TRUE, 
                    indv_points = TRUE, legend = FALSE)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/Mutation_spectrum_cfDNA_postRT_.pdf", width = 12, height = 5)
plot_post<-grid.arrange(p1_post, p2_post, p3_post, ncol = 3)
print(plot_post)
dev.off()


#p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)
#
#p5 <- plot_spectrum(type_occurrences, CT = TRUE, 
#                    legend = TRUE, error_bars = "stdev")
#grid.arrange(p4, p5, ncol = 2, widths = c(4, 2.3))


## 96 muttaional profile
mut_mat <- mut_matrix(vcf_list = sar_grl, ref_genome = ref_genome)
head(mut_mat)


plot_96_profile(mut_mat[, c(1: 7)])
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/pyclone_trunck_vs_branch_Mutation_spectrum_cfDNA_96profile.pdf", width = 10, height = 40)
plot_96<-plot_96_profile(mut_mat[, c(1:38)])
print(plot_96)
dev.off()


#### mutational signature
#### denovo signature
# 1: NMF
mut_mat <- mut_mat + 0.0001
library("NMF")
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/pyclone_trunck_vs_branch_NMF_optimal_factorization_rank_cfDNA_96profile.pdf", width = 5, height = 5)
estimate<-plot(estimate)
plot(estimate)
print(estimate)
dev.off()

### I tool three signature (the most common approach is to choose the smallest rank for which cophenetic correlation coefficient starts decreasing)
nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 10, single_core = TRUE)

nmf_res_contribution<-data.frame(nmf_res$contribution)
write.table(nmf_res_contribution, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/out_res_pyclone_trunck_vs_branch_nmf_res_contribution_cfDNA_96profile.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

#combi_mat = rbind(indel_counts, dbs_counts)
#nmf_res_combi <- extract_signatures(combi_mat, rank = 2, nrun = 10, single_core = TRUE)

#B: Bayesian NMF
# BiocManager::install("ccfindR")

### get similar signatures with Cosmic
signatures = get_known_signatures()
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)

#Visualizing NMF results

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/pyclone_trunck_vs_branch_96profile_noindel_NMF_cfDNA.pdf", width = 7, height = 7)
sig_profile<-plot_96_profile(nmf_res$signatures, condensed = TRUE)
print(sig_profile)
dev.off()


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/pyclone_trunck_vs_branch_barplot_contribution_noindel_NMF_cfDNA.pdf", width = 5, height = 5)
contribution<-plot_contribution(nmf_res$contribution, nmf_res$signature,
  mode = "relative"
)
print(contribution)
dev.off()


### costumized barplot
freq<-type_occurrences_good[,1:6]
freq<-data.frame(rowSums(freq))
freq$sample_id <- row.names(freq)
freq$RT_status<-gsub(".*_","",freq$sample_id)
rownames(freq)<-NULL

denovo_exp<-data.frame(t(nmf_res_contribution))
colnames(denovo_exp)<-c("SBS1-like","SBSA","SBSB")
samples<-rownames(denovo_exp)

our_res_exp<-NULL
for (jj in 1:length(samples)){

     denovo_exp_focal<-data.frame(t(data.frame(denovo_exp[(rownames(denovo_exp)[jj]),])))
     denovo_exp_focal$denovo_sig<-rownames(denovo_exp_focal)
     rownames(denovo_exp_focal)<-NULL
     freq_foc<-freq[freq$sample_id==samples[jj],]
     denovo_exp_focal$relative_contribution<-(denovo_exp_focal[,1]*freq_foc$rowSums.freq.)
     denovo_exp_focal$sample_id<-samples[jj]
     denovo_exp_focal$RT_status<-gsub(".*_","",samples[jj])
     colnames(denovo_exp_focal)<-c("contribution","denovo_sig","relative_contribution","sample_id","RT_status")
     our_res_exp<-rbind(denovo_exp_focal,our_res_exp)
}

rownames(our_res_exp)<-NULL
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/pyclone_trunck_vs_branch_barplot_cexp_contribution_customized_cfDNA.pdf", width = 5, height = 5)
plota<-ggplot(our_res_exp, aes(x = sample_id, y= contribution,fill = denovo_sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   

colnames(our_res_exp)<-c("contribution","denovo_sig","relative_contribution","sample_id","pyclone")

##### compare if there is any difference between signatures among different pre,post and noRT samples ###
SBS1.like<-our_res_exp[our_res_exp$denovo_sig=="SBS1.like",]
SBSA<-our_res_exp[our_res_exp$denovo_sig=="SBSA",]
SBSB<-our_res_exp[our_res_exp$denovo_sig=="SBSB",]



### box plots of difference betwen groups in each signature ###

plot_sig1<-ggplot(SBS1.like, aes(x=pyclone, y=contribution, fill = pyclone)) + 
  geom_boxplot(outlier.colour="white",
                outlier.size=2) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="denovo_exposure_SBS1.like")
wilcox.test(contribution~ pyclone,data = SBS1.like,exact = FALSE) 


plot_sig2<-ggplot(SBSA, aes(x=pyclone, y=contribution, fill = pyclone)) + 
  geom_boxplot(outlier.colour="white",
                outlier.size=2) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
                labs(title="denovo_exposure_SBSA")
wilcox.test(contribution~ pyclone,data = SBSA,exact = FALSE) 


plot_sig3<-ggplot(SBSB, aes(x=pyclone, y=contribution, fill = pyclone)) + 
  geom_boxplot(outlier.colour="white",
                outlier.size=2) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
                labs(title="denovo_exposure_SBSB")
wilcox.test(contribution~ pyclone,data = SBSB,exact = FALSE) 


library("gridExtra")
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/boxplot_pyclone_trunck_branch_cexp_contribution_customized_cfDNA.pdf", width = 4, height = 8)
plot_all_box<-grid.arrange(plot_sig1, plot_sig2,plot_sig3, nrow= 3)
print(plot_all_box)
dev.off()




#########################
### plot the heatmap ####
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/heatmap_pyclone_trunck_branch_contribution_noindel_NMF_cfDNA.pdf", width = 6, height = 6)
heatmap_contribution<-plot_contribution_heatmap(nmf_res$contribution, 
                          cluster_samples = TRUE, 
                          cluster_sigs = TRUE)
print(heatmap_contribution)
dev.off()


hclust_signatures <- cluster_signatures(nmf_res$signatures, method = "average")
signatures_order <- colnames(nmf_res$signatures)[hclust_signatures$order]
signatures_order


hclust_samples <- cluster_signatures(mut_mat, method = "average")
samples_order <- colnames(mut_mat)[hclust_samples$order]
samples_order



plot_contribution_heatmap(nmf_res$contribution,
  sig_order = signatures_order, sample_order = samples_order,
  cluster_sigs = FALSE, cluster_samples = FALSE
)



plot_compare_profiles(mut_mat[, 1],
  nmf_res$reconstructed[, 1],
  profile_names = c("Original", "Reconstructed"),
  condensed = TRUE
)


####################################
#### Cosmic Signature Refitting ####
#### Find mathematically optimal contribution of COSMIC signatures
fit_res <- fit_to_signatures(mut_mat, signatures)
fit_contribution<-data.frame(fit_res$contribution)
write.table(fit_contribution, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/out_res_pyclone_trunck_branch_fit_to_cosmic_contribution_signature_cfDNA.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/barplot_pyclone_trunck_branch_fit_to_cosmic_contribution_signature_cfDNA.pdf", width = 10, height = 10)
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
   samples[ss]
   colnames(fit_contribution_focal)<-c("contribution","cosmic")
   fit_contribution_focal_order<-fit_contribution_focal[order(fit_contribution_focal$contribution,decreasing = TRUE),]   
   fit_contribution_focal_order<-fit_contribution_focal_order[fit_contribution_focal_order$contribution>0,]
   fit_contribution_focal_order<-fit_contribution_focal_order[1:5,]
   fit_contribution_focal_order
   fit_contribution_focal_order$sample_id<-rep(samples[ss])
   out_res_cos<-rbind(out_res_cos,fit_contribution_focal_order)
}

rownames(out_res_cos)<-NULL
out_res_cos<-out_res_cos[out_res_cos$sample_id!="cosmic",]
out_res_cos$contribution<-as.integer(out_res_cos$contribution)


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/barplot_cosmic_contribution_customized_cfDNA.pdf", width = 5, height = 5)
plotcosmic<-ggplot(out_res_cos, aes(x = sample_id, y= contribution,fill = cosmic)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plotcosmic)
dev.off() 


#### unique signatires to trunck vs. branch mutations

trunck_sigs<-out_res_cos[grepl("trunck",out_res_cos$sample_id),]
trunck_sigs_ord<-trunck_sigs[order(trunck_sigs$contribution, decreasing = TRUE),]


branch_sigs<-out_res_cos[grepl("branch",out_res_cos$sample_id),]
branch_sigs_ord<-branch_sigs[order(branch_sigs$contribution, decreasing = TRUE),]


trunck_sigs[trunck_sigs$cosmic%in%branch_sigs$cosmic,]



###pheatmap(out_res_cos,cluster_rows =FALSE)


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
sbs3<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS3",]
sbs8<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS8",]
sbs7a<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS7a",]
sbs15<-out_res_cos2_pre_post[out_res_cos2_pre_post$cosmic=="SBS15",]

#<-out_res_cos_pre_post[out_res_cos_pre_post$cosmic=="SBSB",]


cosmic1<-ggplot(sbs1, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="cosmic SBS1")
wilcox.test(contribution~ RT_status,data = sbs1,exact = FALSE) 
                   
                                  


cosmic13<-ggplot(sbs13, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="cosmic SBS13")
wilcox.test(contribution~ RT_status,data = sbs13,exact = FALSE) 



cosmic3<-ggplot(sbs3, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="cosmic SBS3")
wilcox.test(contribution~ RT_status,data = sbs3,exact = FALSE) 


cosmic8<-ggplot(sbs8, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="cosmic SBS8")
wilcox.test(contribution~ RT_status,data = sbs8,exact = FALSE) 


cosmic7a<-ggplot(sbs7a, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="cosmic SBS7a")
wilcox.test(contribution~ RT_status,data = sbs7a,exact = FALSE) 



cosmic15<-ggplot(sbs15, aes(x=RT_status, y=contribution, fill = RT_status)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=4) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
                labs(title="cosmic SBS15")
wilcox.test(contribution~ RT_status,data = sbs15,exact = FALSE) 


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/cfDNA/boxplot_cosmic_contribution_between_pre_post_compsrison_cfDNA.pdf", width = 8, height = 12)
plot_cosmics<-grid.arrange(cosmic1,cosmic13,cosmic3,cosmic8,cosmic7a,cosmic15,ncol=2)
print(plot_cosmics)
dev.off()

##Stricter refitting
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
cos_sim_samples_

