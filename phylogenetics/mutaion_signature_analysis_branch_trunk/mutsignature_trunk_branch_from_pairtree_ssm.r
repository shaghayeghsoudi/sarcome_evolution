

## analyze mutational signature in branch vs. trunk with mutSignatures from ssm files


#for f in *.tsv; do mv -- "$f" "${f%.tsv}_noRT.tsv"; done
#for f in *.vcf; do (echo "##fileformat=VCFv4.1" && cat "$f") > "${f%.vcf}_corrected.vcf" && rm "$f"; done 
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

###################################
### load metadata ###
################################### 
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")

meta$RTstatus<-ifelse(meta$sequenceofsamplevRT== "beforeRT", "noRT",
                        ifelse(meta$sequenceofsamplevRT== "afterRT", "postRT",
                        ifelse(meta$sequenceofsamplevRT== "nopreopRT", "noRT",
                        "-")))

meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)
meta$identifier<-paste(meta$unique_sample_id,meta$RTstatus, sep = "_")

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


write.table(tree_mutation, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/data/out_res_muttaions_assigned_on_branch_trunk_based_on_pairtree.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

tree_mutation<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/data/out_res_muttaions_assigned_on_branch_trunk_based_on_pairtree.txt", header = TRUE)


tree_mutation$TB_id2<-gsub("trunk_founder","early",tree_mutation$TB_id1)   ## set id you want to run on trunk vs. branch or founder vs nonfounder
tree_mutation$TB_id2<-gsub("branch","late",tree_mutation$TB_id2)


tree_mutation_good<-tree_mutation[,c("Chrom","Pos","Ref","Alt","TB_id1")]
muts<-c("A","T","C","G")
input_sigmutation_snv<-tree_mutation_good[tree_mutation_good$Ref%in%muts,]
input_sigmutation_snv<-input_sigmutation_snv[input_sigmutation_snv$Alt%in%muts,]

colnames(input_sigmutation_snv)<-c("CHROM","POS","REF","ALT","SAMPLEID")

table_count<-table(input_sigmutation_snv$SAMPLEID)
input_sigmutation_snv <- input_sigmutation_snv[input_sigmutation_snv$SAMPLEID %in% names(table_count[table_count >=10]), ]
input_sigmutation_snv$SAMPLEID<-gsub("preRT","noRT",input_sigmutation_snv$SAMPLEID)


input_sigmutation_snv_pre<-input_sigmutation_snv[grepl("noRT",input_sigmutation_snv$SAMPLEID),]  ### make a subset of only no-radiation samples
########################################################################
########################################################################
#### start with Mutsignatures ####
#### extrcat preRT samples
#y_pre<-input_sigmutation_snv[grepl("noRT",input_sigmutation_snv$SAMPLEID),]
y_pre<-input_sigmutation_snv_pre
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


write.table(y_pre, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/trunk_branch/trunk_branch_two_denovo_and_COSMIC_mut10/outres_attachMutType_onlynoRT_trunk_branch_based_on_pairtree_threshold10mut.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

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


#Downstream analyses and visualization
## examine the results
# Retrieve signatures (results)
pre.sig <- pre.analysis$Results$signatures

# Retrieve exposures (results)
pre.exp <- pre.analysis$Results$exposures


# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
for (pp in 1:2){

  pdf(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/trunk_branch/trunk_branch_two_denovo_and_COSMIC_mut10/denovo_MutSignature_onlynoRT_sign_",pp,"_based_on_pairtree_denovotwo.pdf",sep = ""))
  msigPlot(pre.sig, signature = pp, ylim = c(0, 0.10))
  dev.off()
}

#png(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/from_new_metadata_Feb2023/barplot_denovo_MutSignature_counts_per_sample_all_samples_merged_regions_boot500.png")
#msigPlot(pre.exp) + 
#  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))
#dev.off()

# Export Signatures as data.frame
xprt <- coerceObj(x = pre.exp, to = "data.frame") 

#head(xprt) %>% kable() %>% kable_styling(bootstrap_options = "striped")
write.table(xprt, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/trunk_branch/trunk_branch_two_denovo_and_COSMIC_mut10/out_put_mutsignatures_xprt_denove_2_signature_onlynoRT_count_based_on_pairtree.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#xprt<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/MutSignatures/Feb_2023/out_put_mutsignatures_xprt_denove_3_signature_count_All_RT.txt", header = TRUE, sep = "\t")


##########
######### do paired t-test of raw count
xprt_good<-data.frame(t(xprt))
xprt_good$sample_tree<-rownames(xprt_good)
xprt_good_t<-xprt_good %>% separate(sample_tree,c("subject","RT","positio"))
xprt_good_t$id<-paste(xprt_good_t$subject ,xprt_good_t$RT, sep = "_")

find_paires<-table(xprt_good_t$subject)
paired_samples<-xprt_good_t[xprt_good_t$subject%in%names(find_paires[find_paires >=2]),]
xprt_good_paired<-xprt_good_t[xprt_good_t$subject %in%paired_samples$subject,]
sigs<-grep("Sign", names(xprt_good_paired), value=TRUE)


out_res<-NULL
for(ss in 1:length(sigs)){

  xprt_focal<-xprt_good_paired[,c(colnames(xprt_good_paired)==sigs[ss],"positio")]
  res_tpaired <- t.test( Sign.01 ~ positio, data =  xprt_focal, paired = TRUE)
  res_tpaired_df<-data.frame( map_df(list(res), tidy))
  out_res<-rbind(res_tpaired_df,out_res)
}


########################################################
##### costom visualization of denovo exposures #########
### stacked bar plot all samples to gether

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

write.table(our_res_exp, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/founder_branch/noRT_twodenovo/out_put_mutsignatures_only_noRT_xprt_relative_contribution_percentage_denove_2_signature_count_trunk_and_branch_based_on_pairtrer.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/founder_branch/noRT_twodenovo/barplot_trunk_branch_2denovo_sigs_onlynoRT_relative_contribution_based_on_pairtree.pdf", width = 6, height = 7)
plota<-ggplot(our_res_exp, aes(x = sample_id, y= relative_contribution,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   

### make boxplots to comapre relative contribution in trunk versus branch for each denovo mutations
our_res_exp$tree_position<-gsub(".*_","",our_res_exp$sample_id)
denovos<-unique(our_res_exp$sig)

my.plot.denovo <- vector(mode = "list", length = 2)  ### adjust based on the number of Cosmic signatures taken

for(dd in 1:length(denovos)){
  focal<-our_res_exp[our_res_exp$sig==denovos[dd],]

  plot_sig1<-ggplot(focal, aes(x=tree_position, y=relative_contribution, fill = tree_position)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
                labs(title=denovos[dd]) 
                 F1<-plot_sig1+ stat_compare_means(method = "t.test")  ## add P-value
                #F1<-plot_sig1 + stat_compare_means()
                #F1 <- plot_sig1 + stat_compare_means(method = "t.test")
                my.plot.denovo[[dd]] <- F1
}


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/founder_branch/noRT_twodenovo/boxplot_2denovo_relative_contribution_onlynoRT_percent_trunk_branch_based_on_pairtree.pdf", width = 4, height = 8)
plots_denovo_all<-multiplot(plotlist = my.plot.denovo[1:3],cols= 1) 
print(plots_denovo_all)
dev.off()


##### make violin plot #####
our_res_exp %>%

  ggplot(aes(fill=tree_position, y=relative_contribution, x=sig)) + 
    geom_violin()+
    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
    #scale_fill_viridis(discrete=T, name="") +
    #theme_ipsum()  +
    xlab("denovo signature") +
    geom_point(position=position_jitterdodge())+
    ylim(0,100) + stat_compare_means(method = "t.test",label = "p.signif", size = 8)
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/founder_branch/noRT_twodenovo/violinplot_2denovo_relative_contribution_onlynoRT_percent_trunk_branch_based_on_pairtree.pdf", width = 6, height = 6)

### exclude postRT samples 
#our_res_exp_naive<-our_res_exp[grep("noRT|preRT",our_res_exp$sample_id),]
#
#our_res_exp_naive %>%
#
#  ggplot(aes(fill=tree_position, y=relative_contribution, x=sig)) + 
#    geom_violin()+
#    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
#    #scale_fill_viridis(discrete=T, name="") +
#    #theme_ipsum()  +
#    xlab("denovo signature") +
#    geom_point(position=position_jitterdodge())+
#    ylim(0,100) + stat_compare_means(method = "t.test",label = "p.signif", size = 8)
#ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/trunk_branch/trunk_branch_two_denovo/violinplot_treatment_naive_3denovo_relative_contribution_percent_trunk_branch_based_on_pairtree.pdf", width = 6, height = 6)

#############################################
### make stacked bar plot for each sample ####
our_res_exp$subject_id<-sub('^([^_]+_[^_]+).*', '\\1', our_res_exp$sample_id)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/founder_branch/noRT_twodenovo/barplot_stacked_onlynoRT_per_sample_2denovo_relative_contribution_trunk_branch_based_on_ppairtree.pdf", width = 14, height  = 8)
plot_stack<-ggplot(our_res_exp,                         # Draw barplot with grouping & stacking
       aes(x = tree_position,
           y = relative_contribution,
           fill = sig)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  facet_wrap(~ subject_id, nrow = 2)
print(plot_stack)
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
rownames(out_res_final_xprt)<-c("Sign.01","Sign.02")
out_res_final_xprt_df<-data.frame(out_res_final_xprt)

write.table(out_res_final_xprt_df, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/trunk_branch/trunk_branch_two_denovo/out_put_mutsignatures_xprt_relative_contribution_onlynoRT_percentage_denove_2_signature_count_trunk_and_branch_based_on_pairtree.txt", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/trunk_branch/trunk_branch_two_denovo/heatmap_2denovo_onlynoRT_frequencies_trunk_branch_bsed_on_pairtree.pdf", width = 7, height = 7)
pheatmap(out_res_final_xprt_df,cluster_rows =FALSE)
dev.off()

### make ordered heatmap ###
#ordered<-c("SRC125_preRT","SRC127_preRT","SRC167_preRT","SRC169_preRT","SRC170_preRT","SRC171_preRT","TB9051_preRT","SRC125_postRT","SRC127_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB22446_postRT","TB9051_postRT","RC130_noRT","SRC168_noRT" ,"SRC172_noRT","SRC173_noRT", "TB13092_noRT","TB13712_noRT","TB13959_noRT","TB9573_noRT")
#out_res_final_xprt_ordered<-out_res_final_xprt_df[order(match(colnames(out_res_final_xprt_df), ordered))]

#pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_ordered_Unclustered.pdf", width = 5, height = 5)
#pheatmap(out_res_final_xprt_ordered,cluster_cols=FALSE)
#dev.off()


########################################################
####### do paired t-test ############

our_res_exp_sig1<-our_res_exp[our_res_exp$sig=="Sign.02",]
our_res_exp_sig1$uniq_subject<-gsub("_.*$","",our_res_exp_sig1$subject_id)

#our_res_exp_sig1_paired<-our_res_exp_sig1[table(our_res_exp_sig1$uniq_subject)>1,]
paired<-c("SRC125","SRC127","SRC171","SRC172","TB13092" ,"TB13712")
our_res_exp_sig1_paired<-our_res_exp_sig1[our_res_exp_sig1$uniq_subject%in%paired,]
res <- t.test( relative_contribution ~ tree_position, data = our_res_exp_sig1_paired, paired = TRUE)

########################################################
############ chi square test of independence ###########
#### test to see if there is asscociation between the signature contribution and tree location (chi square test)
samples_TB<-c("TB13712_noRT","TB13712_noRT","TB13092_noRT","TB13092_noRT","SRC172_noRT" ,"SRC172_noRT","SRC171_noRT" ,"SRC171_noRT" ,"SRC127_noRT" ,"SRC127_noRT", "SRC125_noRT", "SRC125_noRT")

for(zz in 1:length(samples_TB)){
  focal_clumn<-out_res_final_xprt_df[,(grep(samples_TB[zz],colnames(out_res_final_xprt_df)))]
  chisq.test(focal_clumn)
}


trunk<-out_res_final_xprt_df[,grep("trunk",colnames(out_res_final_xprt_df))]
rowsum_trunk<-rowSums(trunk)
branch<-out_res_final_xprt_df[,grep("branch",colnames(out_res_final_xprt_df))]
rowsum_branch<-rowSums(branch)

rowsum_both<-rbind(rowsum_trunk,rowsum_branch)
rownames(rowsum_both)<-c("trunk","branch")


our_res_denovo_test<-NULL
for(ll in 1:nrow(rowsum_both)){
  focal_row<-t(data.frame(rowsum_both[ll,]))
  rownames(focal_row)<-rownames(rowsum_both)[ll]
  focal_row<-(focal_row/rowSums(focal_row))*100
  our_res_denovo_test<-rbind(focal_row,our_res_denovo_test)
}

test <- chisq.test(our_res_denovo_test)

###############################################################
###############################################################
# Retrieve COSMIC signatures from online repo, and then subset
###############################################################
###############################################################
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
#cosmx<-cosmix[c(1, 2,3, 5,6,13)]
cosmx<-cosmix[c(1,2,5,13)]


# match OVcar and COSMIC signatures
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/trunk_branch/trunk_branch_two_denovo/heatmap_cusine_similarities_onlynoRT_2denovo_cosmic_all_samples.pdf", width = 7, height = 7)
mSign.sar <- matchSignatures(mutSign = pre.sig, reference = cosmix)
print(mSign.sar$plot)
dev.off()


blca.expo1 <- resolveMutSignatures(mutCountData = xx,       #### cosmic sigs
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx,      #### 
                                   signFreqData = pre.sig)

blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.1x_df <- coerceObj(x = blca.exp.1x, to = "data.frame") 
write.table(blca.exp.1x_df, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/out_put_mutsignatures_blca.exp.1x_cosmic_signature_count4_trunk_branch_based_on_pairtree_founder.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


#################
###### do paired t-test of raw cosmic counts
blca.exp.1x_good<-data.frame(t(blca.exp.1x_df))
blca.exp.1x_good$sample_tree<-rownames(blca.exp.1x_good)
blca.exp.1x_good_t<-blca.exp.1x_good %>% separate(sample_tree,c("subject","RT","positio"))
blca.exp.1x_good_t$id<-paste(blca.exp.1x_good_t$subject ,blca.exp.1x_good_t$RT, sep = "_")
#paired<-c("SRC125","SRC127","SRC171","SRC172","TB13092" ,"TB13712")
paired<-c("SRC127","SRC169","SRC170","SRC171","TB9051")
blca.exp.1x_good_t_paired<-blca.exp.1x_good_t[blca.exp.1x_good_t$subject %in%paired,]
blca.exp.1x_good_t_paired_sig1<-blca.exp.1x_good_t_paired[,c("COSMIC.13","positio")]
res <- t.test( COSMIC.13 ~ positio, data = blca.exp.1x_good_t_paired_sig1, paired = TRUE)
res






blca.exp.1x.freq <- blca.expo1$Results$freq.result
blca.exp.1x_df.freq <- coerceObj(x = blca.exp.1x.freq, to = "data.frame") 
write.table(blca.exp.1x_df.freq, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/out_put_mutsignatures_blca.exp.1x_cosmic_signature_freq4_trunk_branch_based_on_pairtree_founder.txt", col.names = TRUE, row.names = TRUE, sep = "\t",quote = FALSE)


blca.exp.2x <- blca.expo2$Results$count.result
blca.exp.2x_df <- coerceObj(x = blca.exp.2x, to = "data.frame") 
write.table(blca.exp.2x_df, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/out_put_mutsignatures_blca.exp.2x_denove_signature_count_All_RT_based_on_pairtree_founder.txt")



# Plot exposures
cosmic_exp<-data.frame(t(blca.exp.1x_df))

our_res_exp_cosmic<-NULL
for(kk in 1:nrow(cosmic_exp)){

  focal_sample<-cosmic_exp[kk,]
  focal_sample<-(focal_sample/rowSums(focal_sample))*100
  count_table<-data.frame(t(focal_sample))
  count_table$sig<-rownames(count_table)
  count_table$sample_id<-rownames(cosmic_exp[kk,])
  colnames(count_table)<-c("relative_contribution","sig","sample_id")
  rownames(count_table)<-NULL
  our_res_exp_cosmic<-rbind(count_table,our_res_exp_cosmic)
}

write.table(our_res_exp, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/out_put_mutsignatures_xprt_relative_contribution_percentage_cosmic_4signature_count_All_RT.txt",col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)


pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/founder_branch/noRT_twodenovo/barplot_4cosmic_cexp_relative_contribution_percentage_trunk_branch_based_on_pairtree_founder.pdf", width = 6, height = 7)
plota<-ggplot(our_res_exp_cosmic, aes(x = sample_id, y= relative_contribution,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   



#############################
##### make violin plot ######
#############################

our_res_exp_cosmic$tree_position<-gsub(".*_","",our_res_exp_cosmic$sample_id)

our_res_exp_cosmic %>%

  ggplot(aes(fill=tree_position, y=relative_contribution, x=sig)) + 
    geom_violin()+
    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
    #scale_fill_viridis(discrete=T, name="") +
    #theme_ipsum()  +
    xlab("denovo signature") +
    geom_point(position=position_jitterdodge())+
    ylim(0,100) + stat_compare_means(method = "t.test",label = "p.signif", size = 8)
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/violinplot_cosmic_relative_contribution_percent_trunk_branch_based_on_pairtree.pdf", width = 6, height = 6)



### exclude postRT samples
our_res_exp_cosmic_naive<-our_res_exp_cosmic[grep("noRT|preRT",our_res_exp_cosmic$sample_id),]

our_res_exp_cosmic_naive %>%

  ggplot(aes(fill=tree_position, y=relative_contribution, x=sig)) + 
    geom_violin()+
    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
    #scale_fill_viridis(discrete=T, name="") +
    #theme_ipsum()  +
    xlab("denovo signature") +
    geom_point(position=position_jitterdodge())+
    ylim(0,100) + stat_compare_means(method = "t.test",label = "p.signif", size = 8)
ggsave(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/violinplot_treatment_cosmic_relative_contribution_percent_trunk_branch_based_on_pairtree.pdf", width = 6, height = 6)


##### make  stacked barplot ########
our_res_exp_cosmic$subject_id<-sub('^([^_]+_[^_]+).*', '\\1', our_res_exp_cosmic$sample_id)


#nobranch<-c("SRC130","SRC150", "SRC168","TB11985")
#our_res_exp_cosmic_TBboth<-our_res_exp_cosmic[!(our_res_exp_cosmic$subject_id%in%nobranch),]
#our_res_exp_cosmic_TBboth$tree_position<-gsub(".*_","",our_res_exp_cosmic_TBboth$sample_id)


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/founder_branch/noRT_twodenovo/barplot_stacked_per_sample_4Cosmic_relative_contribution_trunk_branch_from_pairtree.pdf", width = 10, height  = 8.5)
plot_stack<-ggplot(our_res_exp_cosmic,                         # Draw barplot with grouping & stacking
       aes(x = tree_position,
           y = relative_contribution,
           fill = sig)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  facet_wrap(~ subject_id, nrow = 3)
  #facet_grid(~ subject_id, nrow = 2)
print(plot_stack)
dev.off()


### boxplots to comapre each signature in branch vs. trunk
#our_res_exp$RTstatus<-gsub("naive","noRT",our_res_exp$RTstatus)
#our_res_exp$RTstatus<-gsub("RTtreatment","postRT",our_res_exp$RTstatus)
names(our_res_exp_cosmic)[1]<-"percent_relative_contribution"

sigs<-unique(our_res_exp_cosmic$sig)
my.plot <- vector(mode = "list", length = 4)  ### adjust based on the number of Cosmic signatures taken

for(ss in 1:length(sigs)){
  focal<-our_res_exp_cosmic[our_res_exp_cosmic$sig==sigs[ss],]

  plot_sig1<-ggplot(focal, aes(x=tree_position, y=percent_relative_contribution, fill = tree_position)) + 
  geom_boxplot(outlier.colour="black",
                outlier.size=0.7) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
                labs(title=sigs[ss]) 
                #F1<-plot_sig1 + stat_compare_means()
                F1 <- plot_sig1 + stat_compare_means(method = "t.test")
                my.plot[[ss]] <- F1
}

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/boxplot_cexp_counts_all_samples_sample_annotated_4cosmic_relative_contribution_percentage_denoco_based_on_pairtree_founder.pdf", width = 9, height = 10)
plots_cosmic_all<-multiplot(plotlist = my.plot[1:4],cols= 2) 
print(plots_cosmic_all)
dev.off()


#### heatmap of denovo exposure frequencies ####
#freq<-data.frame(table(y_pre$SAMPLEID))
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

pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/signature_analysis/truck_branch_signature_updated_solution/mutsignature/heatmap_4Cosmic_signatures_frequencies_all_samples_annotated_clustered_based_on_pairtree_founder.pdf", width = 8, height = 8)
plot_heat<-pheatmap(out_res_final_xprtcos_df,cluster_rows =FALSE)
print(plot_heat)
dev.off()

#ordered<-c("SRC125_preRT","SRC127_preRT","SRC167_preRT","SRC169_preRT","SRC170_preRT","SRC171_preRT","TB9051_preRT","SRC125_postRT","SRC127_postRT","SRC167_postRT","SRC169_postRT","SRC170_postRT","SRC171_postRT","TB22446_postRT","TB9051_postRT","RC130_noRT","SRC168_noRT" ,"SRC172_noRT","SRC173_noRT", "TB13092_noRT","TB13712_noRT","TB13959_noRT","TB9573_noRT")
#out_res_final_xprt_ordered<-out_res_final_xprt_df[order(match(colnames(out_res_final_xprt_df), ordered))]

#pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/MutationalSignature/cfDNA/heatmap_denovo_frequencies_cfDNA_sample_annotated_ordered_Unclustered.pdf", width = 5, height = 5)
#pheatmap(out_res_final_xprt_ordered,cluster_cols=FALSE)
#dev.off()


trunk_cosmic<-out_res_final_xprtcos_df[,grep("trunk",colnames(out_res_final_xprtcos_df))]
trunk_cosmic_row<-rowSums(trunk_cosmic)
branch_cosmic<-out_res_final_xprtcos_df[,grep("branch",colnames(out_res_final_xprtcos_df))]
branch_cosmic_row<-rowSums(branch_cosmic)

rowsum_both_cosmic<-rbind(trunk_cosmic_row,branch_cosmic_row)
rownames(rowsum_both_cosmic)<-c("trunk","branch")


our_res_cosmic_test<-NULL
for(ll in 1:nrow(rowsum_both_cosmic)){
  focal_row<-t(data.frame(rowsum_both_cosmic[ll,]))
  rownames(focal_row)<-rownames(rowsum_both_cosmic)[ll]
  focal_row<-(focal_row/rowSums(focal_row))*100
  our_res_cosmic_test<-rbind(focal_row,our_res_cosmic_test)
}

test <- chisq.test(our_res_cosmic_test)