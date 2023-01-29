
### Mutational signature analysis with Muttaionalpatterns package in R
## for vignette see: https://bioconductor.org/packages/release/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.html

#rm(list = ls())
library(MutationalPatterns)
library(BSgenome)
library(dplyr)
library(stringr)
library(tidyr)
#head(available.genomes())

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")


meta<-read.table(file = "metadata_updated_shaghayegh_Jan2023.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)

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
mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")

mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

mutaion_RT$RTtretment<-ifelse(mutaion_RT$RT_status== "preRT", "naive",
                        ifelse(mutaion_RT$RT_status== "noRT", "naive",
                        ifelse(mutaion_RT$RT_status== "postRT", "irradiated",
                        "-")))

mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id.x,mutaion_RT$RT_status, sep = "_")
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
#                                 ref_genome, type = "indel")

######
##indel_grl <- get_indel_context(indel_grl, ref_genome)
#head(indel_grl[[1]], n = 3)

#Next count the number of indels per type. This results in a matrix that is similar to the mut_mat matrix.

#indel_counts <- count_indel_contexts(indel_grl)
#head(indel_counts)
#plot the indel spectra
#plot_indel_contexts(indel_counts, condensed = TRUE)

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
write.table(type_occurrences, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/type_occurrences_muutaionpattern_exome_noindels_Jan2023.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


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
ggsave(filename="~/Desktop/boxplot_nucleotide_types_per_radiation.jpg", plot=box_type, height = 5, width = 8)


## 96 muttaional profile
mut_mat <- mut_matrix(vcf_list = sar_grl, ref_genome = ref_genome)
head(mut_mat)
plot_96_profile(mut_mat[, c(1: 7)])


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/Mutation_spectrum_exome_96profile_noindel_Jan2023.pdf", width = 10, height = 30)
plot_A<-plot_96_profile(mut_mat[, c(1:27)])
print(plot_A)
dev.off()


#### mutational signature
#### denovo signature
# 1: NMF
mut_mat <- mut_mat + 0.0001
library("NMF")
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/NMF_optimal_factorization_rank_exome_96profile_noindel_Jan2023.pdf", width = 5, height = 5)
estimate<-plot(estimate)
plot(estimate)
print(estimate)
dev.off()

### I tool three signature (the most common approach is to choose the smallest rank for which cophenetic correlation coefficient starts decreasing)
nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 10, single_core = TRUE)

nmf_res_contribution<-data.frame(nmf_res$contribution)
rownames(nmf_res_contribution)<-c("SBS1-like","SBSA","SBSB")
write.table(nmf_res_contribution, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/out_res_nmf_res_contribution_exome_96profile_noindel_Jan2023.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

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


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/heatmap_contribution_noindel_NMF_exome.pdf", width = 7, height = 5)
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


######
#Signature refitting
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
cos_sim_samples_signatures[1:3, 1:3]                    