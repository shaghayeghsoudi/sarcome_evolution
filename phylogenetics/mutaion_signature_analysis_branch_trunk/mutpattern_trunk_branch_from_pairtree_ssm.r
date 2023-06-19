## analyze mutational signature in branch vs. trunk with mutSignatures


#for f in *.tsv; do mv -- "$f" "${f%.tsv}_noRT.tsv"; done
#for f in *.vcf; do (echo "##fileformat=VCFv4.1" && cat "$f") > "${f%.vcf}_corrected.vcf" && rm "$f"; done 


rm(list = ls())
library(MutationalPatterns)
library(BSgenome)
library("rjson")
library(ggplot2)
library(stringr)
library(ggpubr)
library(pheatmap)
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

source("http://peterhaschke.com/Code/multiplot.R")

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

###################################
### load metadata ###
################################### 
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")

meta$RTstatus<-ifelse(meta$sequenceofsamplevRT== "beforeRT", "noRT",
                        ifelse(meta$sequenceofsamplevRT== "afterRT", "postRT",
                        ifelse(meta$sequenceofsamplevRT== "nopreopRT", "noRT",
                        "-")))

meta$sequenceofsamplevRT<-gsub("beforeRT","preRT",meta$sequenceofsamplevRT)
meta$sequenceofsamplevRT<-gsub("afterRT","postRT",meta$sequenceofsamplevRT)
meta$sequenceofsamplevRT<-gsub("nopreopRT","noRT",meta$sequenceofsamplevRT)

meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)
#mutaion_RT<-merge(py_data_good,meta, by.x = "sample_id", by.y = "sampleid")

meta$identifier<-paste(meta$unique_sample_id,meta$RTstatus, sep = "_")
#mutaion_RT$pos_identifier<-gsub("chr","",mutaion_RT$pos_id)
#samples<-unique(mutaion_RT$unique_sample_id.x)

###################################
### load tree Json files ###
###################################
## read tree location file
tree<-read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/paired_cluster_location_from_pirtree.csv", header = TRUE, sep = ",")
#tree<-tree[tree$sample_id!="TB12052_postRT",]
tree$sample_id <- gsub("SRC168_preRT","SRC168_noRT",tree$sample_id)
samples<-unique(tree$sample_id)

tree<-tree[!(tree$cluster_id==0),]
tree$sample_cluster<-paste(tree$sample_id ,tree$cluster_id , sep = "_")


############################################################
############################################################
################# assign muttaions to clusters ############
### load pyc.json files 
file_pyc<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/April_updated1/full/pyc.jsons", pattern = "*pyc.out.json", full.names = TRUE)

out_res_pp<-NULL
for (pp in 1:length(file_pyc)){
     
     focal_file_pyc<-file_pyc[pp] 
     json_focal_pyc <- fromJSON(paste(readLines(focal_file_pyc), collapse=""))  ## parse Json tree file into R
     ## adjust samples name to match with variant file
     subject_id<-gsub(".pyc.out.json", "",sub(".*/", "", focal_file_pyc)) 
     subject_id<-gsub("-","_",subject_id)
     
     
     cluster<-json_focal_pyc$clusters
     out_res_jj<-NULL
      for(jj in 1:length(cluster)){

          cluster_focal<-data.frame("id"=cluster[jj])
          cluster_focal$cluster_id<-jj
          cluster_focal$subject_id<-subject_id
          names(cluster_focal)[1]<-"id"
          out_res_jj<-rbind(cluster_focal,out_res_jj)
      }
out_res_pp<-rbind(out_res_jj,out_res_pp)

}
     
out_res_pp$checkpoint<-paste(out_res_pp$subject_id ,out_res_pp$id , sep = "_")    
#    id cluster_id   subject_id
#1  s2          8 TB9051_preRT
#2  s7          8 TB9051_preRT
#3 s10          8 TB9051_preRT
#4 s14          8 TB9051_preRT    
  
#focal_file_tree<-file_trees[tt] 
#son_focal <- fromJSON(paste(readLines(focal_file_tree), collapse=""))  ## parse Json tree file into R



############################################################
############################################################
### load  files to make pyclone muttaion like file ######
############################################################
############################################################
filenames_py <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/April_updated1/full/ssms",pattern="*.ssm", full.names = TRUE)  ### load pyclone cluster files
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")[,c("name","id","mutation_id")]
     })

for (i in 1:length(attackStats_py)){
    attackStats_py[[i]]<-cbind(attackStats_py[[i]],filenames_py[i])
    colnames(attackStats_py[[i]])<-c("position","id","mutation_id","file_name")
    }

py_data <- do.call("rbind", attackStats_py) 
py_data$subject_id<-sub('.*\\/', '', py_data$file_name)
py_data$subject_id<-gsub(".ssm","",py_data$subject_id)
py_data$subject_id<-gsub("-","_",py_data$subject_id)   # "SRC125_postRT"

py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))  ## make seperate columns for chrom, pos, ref and alt allels
#py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
#py_data_good$subject_id<-gsub("_.*$","",py_data_good$sample_id)

py_data_good$position<-gsub("chr","",py_data_good$position)
colnames(py_data_good)[colnames(py_data_good) == "position"] <- "pos_id"
py_data_good<-py_data_good[,-7]

py_data_good$checkpoint<-paste(py_data_good$subject_id ,py_data_good$id , sep = "_")    


mut_table<-merge(out_res_pp,py_data_good,by.x = "checkpoint", by.y= "checkpoint")[,c("Chrom","Pos","Ref","Alt","subject_id.y","checkpoint","id.x","cluster_id")]
mut_table$sample_cluster<-paste(mut_table$subject_id.y ,mut_table$cluster_id ,sep = "_")

tree_mutation<-merge(mut_table,tree, by.x = "sample_cluster",by.y = "sample_cluster")

tree_mutation$TB_id1<-paste(tree_mutation$sample_id,tree_mutation$tree_location,sep = "_")
tree_mutation$TB_id2<-paste(tree_mutation$sample_id,tree_mutation$tree_detail, sep = "_")


#out_res_good<-out_res[,c("Chrom","Pos","Ref","Alt","TB_id")]
#colnames(out_res_good)<-c("CHROM","POS","REF","ALT","SAMPLEID")

tree_mutation$ID<-rep(".",nrow(tree_mutation))
tree_mutation$QUAL<-rep(".",nrow(tree_mutation))
tree_mutation$FILTER<-rep(".",nrow(tree_mutation))
tree_mutation$FORMAT<-rep("AD:DP",nrow(tree_mutation))
#tree_mutation$Ref_count<-(tree_mutation$TumorDepth-tree_mutation$TumorReads)
tree_mutation$Ref_count<-0
tree_mutation$allele<-paste(".",(paste(tree_mutation$TumorReads,tree_mutation$Ref_count, sep =":")),sep=",")
tree_mutation$Ref_count<-paste("t_ref_count=",tree_mutation$Ref_count,sep = "")
tree_mutation$TumorReads<-paste("t_alt_count=",tree_mutation$TumorReads,sep = "")
tree_mutation$INFO<-paste(tree_mutation$TumorReads,tree_mutation$Ref_count, sep =";" )

#out_res$Chromosome<-gsub("chr","",out_res$Chromosome)
#muts<-c("A","T","C","G")
#out_res<-out_res[out_res$Ref%in%muts,]
#out_res<-out_res[out_res$Alt%in%muts,]

#out_res$Chromosome<-gsub("chr", "",out_res$Chromosome)
tree_mutation$id1<-paste(tree_mutation$Chrom,tree_mutation$Pos, sep = ":")
tree_mutation$id2<-paste(tree_mutation$Ref,tree_mutation$Alt, sep = "/")
tree_mutation$ID<-paste(tree_mutation$id1,tree_mutation$id2, sep = "_")



muttaion_selected_cols<-tree_mutation[,c("Chrom","Pos","ID","Ref","Alt","QUAL","FILTER","INFO","FORMAT","allele","TB_id2")]
colnames(muttaion_selected_cols)<-c("#CHROM","POS","ID","REF","ALT"	,"QUAL"	,"FILTER",	"INFO",	"FORMAT","SAMPLE", "sample_id")
colnames(muttaion_selected_cols)<-c("#CHROM","POS","ID","REF","ALT"	,"QUAL"	,"FILTER",	"INFO",	"FORMAT","SAMPLE", "sample_id")

samples<-unique(muttaion_selected_cols$sample_id)


### make and save vcf format files ###
for (ii in 1:length(samples)){

    focal_vcf<-muttaion_selected_cols[muttaion_selected_cols$sample_id==samples[ii],]
    colnames(focal_vcf)[10]<-samples[ii]
    focal_vcf<-focal_vcf[,-11]
    #focal_vcf<- focal_vcf[!grepl("-",focal_vcf$REF),]
    #focal_vcf<- focal_vcf[!grepl("-",focal_vcf$ALT),]
    write.table(focal_vcf, file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_vcf/",samples[ii],".filtered.variants.foundertrunk_branch.May2023.vcf", sep = ""),row.names= FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

##############################
##### run mutpatterns ########
vcf_files <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_vcf/unradiatted", pattern = "*.filtered.variants.foundertrunk_branch.May2023.vcf",full.names = TRUE)
#sample_names <- c("SRC125_postRT","SRC125_noRT","SRC127_postRT","SRC127_noRT","SRC130_noRT","SRC150_postRT","SRC167_postRT","SRC167_noRT","SRC168_noRT","SRC169_postRT","SRC169_noRT","SRC170_postRT","SRC170_noRT","SRC171_postRT","SRC171_noRT","SRC172_noRT","SRC173_noRT","TB11985_postRT","TB12052_postRT","TB13092_noRT","TB13712_noRT","TB13959_noRT","TB22446_postRT","TB8016_noRT","TB9051_postRT","TB9051_noRT","TB9573_noRT")
sample_names<-gsub("/Users/shsoudi/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_vcf/","",gsub(".filtered.variants.foundertrunk_branch.May2023.vcf","",vcf_files))

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
write.table(type_occurrences, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/type_occurrences_muutaionpattern_foundertrunk_branch_filt20_May2023.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


############################################################################################
############################################################################################
#### customized plots 
#### make a bar plot of relative contribution of each base substitution in each sample
#type_occurrences$RTstatus<-gsub(".*_","",rownames(type_occurrences))
type_occurrences$time <- sub('.*_\\s*', '', rownames(type_occurrences))

type_occurrences$sample_id<-rownames(type_occurrences)
rownames(type_occurrences)<-NULL

type_occurrences<-type_occurrences[,c("C>A","C>G","C>T","T>A","T>C","T>G","sample_id","time")]

out_res_full<-NULL
for(pp in 1:nrow(type_occurrences)){
  
  focal_pp<-type_occurrences[pp,]
  focal_pp$mut_count<-data.frame(total_muts = apply(focal_pp[1:6], 1, sum))
  focal_pp_propotion<-data.frame(lapply(focal_pp[1:6], function(x) x /as.integer(focal_pp$mut_count)))
  names(focal_pp_propotion)<-paste("prop",names(focal_pp_propotion), sep = "_")
  #focal_pp_propotion$sample_id<-focal_pp$sample_id
  #focal_pp_propotion$RTstatus<-focal_pp$RTstatus
  #focal_count_proportion<-cbind(focal_pp,focal_pp_propotion)
  #types_prop<-grep("prop",names(focal_count_proportion), value= TRUE)
  types_prop<-c("prop_C.A","prop_C.G","prop_C.T","prop_T.A","prop_T.C","prop_T.G")

  my_type<-data.frame("relative_contribution"=t(focal_pp_propotion[,names(focal_pp_propotion)%in%types_prop]))
  my_type$bases<-rownames(my_type)
  rownames(my_type)<-NULL
  my_type$sample_id<-focal_pp$sample_id
  my_type$time<-focal_pp$time
  my_type$bases<-gsub("prop_","",my_type$bases)


  out_res_full<-rbind(my_type,out_res_full)
}

write.table(type_occurrences, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/type_occurrences_relative_contribution_muutaionpattern_trunkfounder_branch_May2023.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)

plot_full_bases<-ggplot(out_res_full, aes(fill=bases , y=relative_contribution, x=sample_id)) + 
    geom_bar(position="fill", stat="identity")
     plot_full_bases + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(filename="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/barplot_single_base_substitution_relative_contribution_foundertrunk_branch_May2023.pdf", height = 7, width = 7)



#### violin plot all full
out_res_full %>%

  ggplot(aes(fill=time, y=relative_contribution, x=bases)) + 
    geom_violin()+
    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
    #scale_fill_viridis(discrete=T, name="") +
    #theme_ipsum()  +
    xlab("denovo signature") +
    geom_point(position=position_jitterdodge())+
    ylim(0,1) + stat_compare_means(method = "t.test",label = "p.signif", size = 8)
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/violinplot_SBS_relative_contribution_trunk_branch_based_on_pairtree.pdf", width = 6, height = 6)


#### violin plot all full but only unradiatted samples
out_res_noRT<-out_res_full[grepl("noRT|preRT",out_res_full$sample_id),]
out_res_noRT %>%

  ggplot(aes(fill=time, y=relative_contribution, x=bases)) + 
    geom_violin()+
    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
    #scale_fill_viridis(discrete=T, name="") +
    #theme_ipsum()  +
    xlab("denovo signature") +
    geom_point(position=position_jitterdodge())+
    ylim(0,1) + stat_compare_means(method = "t.test",label = "p.signif", size = 8)
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/trunk_branch_analysis/violinplot_SBS_relative_contribution_trunk_branch_based_on_pairtree_only_Unradiated.pdf", width = 6, height = 6)



### make a box plot of relative contribution of each base substitution ###
p <- ggplot(data = out_res_full, aes(x=bases, y=relative_contribution)) + 
             geom_boxplot(aes(fill=time)) +
             geom_jitter() 
             p + stat_compare_means(aes(group = time), label = "p.signif")
#p + stat_compare_means(aes(group = time), label = "p.format")
#p + facet_wrap( ~ bases, scales="free")
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/boxplot_SBS_relative_contribution_trunk_branch_based_on_pairtree_full.pdf", width = 6, height = 6)



### make a box plot of relative contribution of each base substitution (only unradiated) ###
p <- ggplot(data = out_res_noRT, aes(x=bases, y=relative_contribution)) + 
             geom_boxplot(aes(fill=time)) +
             geom_jitter() 
             p + stat_compare_means(aes(group = time), label = "p.signif")
#p + stat_compare_means(aes(group = time), label = "p.format")
#p + facet_wrap( ~ bases, scales="free")
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/trunk_branch_analysis/boxplots_SBS_relative_contribution_trunk_branch_based_on_pairtree_only_Unradiated.pdf", width = 6, height = 6)

############################################################################################
############################################################################################
## 96 muttaional profile
mut_mat <- mut_matrix(vcf_list = sar_grl, ref_genome = ref_genome)
head(mut_mat)
plot_96_profile(mut_mat[, c(1: 7)])


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/Mutation_spectrum_trunk_branch_96profile_May2023.pdf", width = 7, height = 35)
plot_A<-plot_96_profile(mut_mat[, c(1:27)])
print(plot_A)
dev.off()

############################################################################################
############################################################################################
#### De novo mutational signature extraction using NMF #
#### denovo signature
# 1: using nonnegative matrix factorization (NMF)
mut_mat <- mut_mat + 0.0001
library("NMF")
estimate <- nmf(mut_mat, rank = 2:6, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/NMF_optimal_factorization_rank_trunk_branch.pdf", width = 5, height = 5)
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
write.table(nmf_res_contribution, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/out_res_nmf_res_contribution_foundertrunk_branch_may2023.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)

png(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/heatmap_contribution_NMF_trunk_branch_full_unclustered.png", width = 700, height = 800)
heatmap_contribution<-plot_contribution_heatmap(nmf_res$contribution,
  #sig_order = signatures_order, sample_order = samples_order,
  cluster_sigs = FALSE, cluster_samples = FALSE
)
print(heatmap_contribution)
dev.off()

### make a customized heatmap of relative contribution of each denovo signature

nmf_res_contribution_prop<-data.frame(lapply(nmf_res_contribution, function(x) x*(1/sum(x))))
rownames(nmf_res_contribution_prop)<-c("SBSA","SBSB","SBS1-like")

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/trunk_branch_analysis/heatmap_contribution_customozed_denovo_founder_notfounder_analysis.pdf", width = 7, height = 7)
custom_heat<-pheatmap((nmf_res_contribution_prop))
print(custom_heat)
dev.off()

#########################################################################
### make a bar plot of relative contribution of each denovo signature ###
out_res_nmf<-NULL
for (jj in 1:ncol(nmf_res_contribution)){

  nmf_focal<-data.frame("contribution"= nmf_res_contribution[,jj])
  normalized<-1/colSums(nmf_focal)
  nmf_focal$relative_contribution<-(nmf_focal$contribution)*normalized
  nmf_focal$signature<-rownames(nmf_res_contribution)
  nmf_focal$patient_ID<-colnames(nmf_res_contribution[jj])
  out_res_nmf<-rbind( nmf_focal,out_res_nmf)

}
write.table(out_res_nmf, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/out_res_normalized_nmf_res_relative_contribution_trunk_branch_May2023.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/96profile_noindel_trunk_branch.pdf", width = 7, height = 7)
plot_96profile<-plot_96_profile(nmf_res$signatures, condensed = TRUE)
print(plot_96profile)
dev.off()


##### make bar plot of relative contribution of each signature
#out_res_nmf$RTstatus<-gsub(".*_","",out_res_nmf$patient_ID)
out_res_nmf$time <- sub('.*_\\s*', '', out_res_nmf$patient_ID)

denovo_contribution<-ggplot(out_res_nmf, aes(fill=signature , y=relative_contribution, x=patient_ID)) + 
    geom_bar(position="fill", stat="identity")
     denovo_contribution + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(filename="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/barplot_mutpatterns_relative_contribution_NMF_denovo_trunkfounder_branch_Full.pdf", height = 7, width = 7)


p <- ggplot(data = out_res_nmf, aes(x=signature, y=relative_contribution)) + 
             geom_boxplot(aes(fill=time)) +
             geom_jitter() 
             p + stat_compare_means(aes(group = time), label = "p.signif")
#p + stat_compare_means(aes(group = time), label = "p.format")
#p + facet_wrap( ~ bases, scales="free")
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/boxplots_SBS_relative_contribution_NMF_denovo_trunkfounder_branch_Full.pdf", width = 6, height = 6)


out_res_nmf %>%

  ggplot(aes(fill=time, y=relative_contribution, x=signature)) + 
    geom_violin()+
    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
    #scale_fill_viridis(discrete=T, name="") +
    #theme_ipsum()  +
    xlab("denovo signature") +
    geom_point(position=position_jitterdodge())+
    ylim(0,1) + stat_compare_means(method = "t.test",label = "p.signif", size = 8)
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/trunk_branch_analysis/violinplot_SBS_relative_contribution_trunk_branch_based_on_pairtree.pdf", width = 6, height = 6)


#############################################
### make stacked bar plot for each sample ####
out_res_nmf$subject_id<-sub('^([^_]+_[^_]+).*', '\\1', out_res_nmf$patient_ID)

#nobranch<-c("SRC130","SRC150", "SRC168","TB11985")
#our_res_exp_TBboth<-our_res_exp[!(our_res_exp$subject_id%in%nobranch),]
#our_res_exp$tree_position<-gsub(".*_","",our_res_exp$sample_id)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/founder_notfounder_analysis/barplot_stacked_per_sample_SBS_relative_contribution_trunk_branch_based_on_ppairtree.pdf", width = 14, height  = 8)
plot_stack<-ggplot(out_res_nmf,                         # Draw barplot with grouping & stacking
       aes(x = time,
           y = relative_contribution,
           fill = signature)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  facet_wrap(~ subject_id, nrow = 2)
print(plot_stack)
dev.off()


#sigs<-c("SBS1-like","SBSA","SBSB","SBSC")
#my.plot <- vector(mode = "list", length = 3)  ### adjust based on the number of Cosmic signatures taken
#
#for(ss in 1:length(sigs)){
#  focal<-out_res_nmf[out_res_nmf$signature==sigs[ss],]
#
#  plot_sig1<-ggplot(focal, aes(x=time, y=relative_contribution, fill = time)) + 
#  geom_boxplot(outlier.colour="black",
#                outlier.size=0.7) +
#                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
#                labs(title=sigs[ss]) 
#                F1<-plot_sig1 + stat_compare_means()
#                #F1 <- plot_sig1 + stat_compare_means(method = "t.test")
#                my.plot[[ss]] <- F1
#}
#
#pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/trunk_branch_analysis/boxplot_relative_contribution_customized_trunk_branch_may2023.pdf", width = 5, height = 12)
#plots_cosmic_all<-multiplot(plotlist = my.plot[1:3],cols= 1) 
#print(plots_cosmic_all)
#dev.off()


### check normality
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/mutational_patterns/summary_all_WES_noindels_regions_added_Feb2023/histogram_denovo_relative_contribution.pdf", width = 7, height = 5)
par(mfrow=c(1,3)) 
hist(SBS1.like$relative_contribution, breaks = 10, main = "SBS1.lik")
hist(SBSA$relative_contribution, breaks = 10, main = "SBSA")
hist(SBSB$relative_contribution, breaks = 10, main = "SBSB")
dev.off()
shapiro.test(SBSB$relative_contribution)


#########################################
#########################################
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

### plot heatmap


nmf_res$contribution


### plot the heatmap (tool made)
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
write.table(fit_contribution, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutpattern/trunk_branch_analysis/out_res_fit_to_cosmic_contribution_May2023.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


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