rm(list = ls())
## load json files 
#aa<-fromJSON(file ="SRC125.tree.json")
#### pyclone outputs ####
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
source("http://peterhaschke.com/Code/multiplot.R")

###
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/trees_from_updated_solutions/data_and_trees")
filenames_py <- list.files("00-inputs-pyclone",pattern="*.pyclone.results.tsv", full.names = TRUE)
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

py_data <- do.call("rbind", attackStats_py) 
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))

py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")


## load meta data 
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")

meta$RTstatus<-ifelse(meta$sequenceofsamplevRT== "beforeRT", "noRT",
                        ifelse(meta$sequenceofsamplevRT== "nopreopRT", "noRT",
                        ifelse(meta$sequenceofsamplevRT== "afterRT", "postRT",
                        "-")))

meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)
paired_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167")
meta_paired<-meta[meta$unique_sample_id%in%paired_samples,]

py_samples<-paired_samples

out_res_py<-NULL
for (zz in 1:length(paired_samples)){

       py_focal<-py_data_good[py_data_good$unique_sample_id==py_samples[zz],]
       ## remove a cluster with only 1 mutation and rename the clusters (new column "new_cluster_id" will be generated)##
       count_table<-table(py_focal$cluster_id)
       to_keep<-data.frame(which(count_table > length(unique(py_focal$sample_id))))
       py_focal_tokeep<-py_focal[py_focal$cluster_id%in%rownames(to_keep),]

       py_focal_tokeep_renamed<-py_focal_tokeep |>
       arrange(cluster_id) |>
       mutate(new_cluster_id = cur_group_id()-1,
         .by = "cluster_id",
         .after = "cluster_id")

       mutaion_RT_focal<-merge(py_focal_tokeep_renamed,meta_paired, by.x = "sample_id", by.y = "sampleid")
       mutaion_RT_focal$identifier<-paste(mutaion_RT_focal$unique_sample_id.x,mutaion_RT_focal$RTstatus, sep = "_")
       mutaion_RT_focal$pos_identifier<-gsub("chr","",mutaion_RT_focal$pos_id)

       ### calulate mean of CP for preRT samples
       focal_mut_pre<-mutaion_RT_focal[grepl("no",mutaion_RT_focal$identifier),]
       focal_mut_pre_mean_cp<-aggregate( cellular_prevalence ~ new_cluster_id, focal_mut_pre, mean )
       focal_mut_pre_mean_cp$RTstatus<-unique(focal_mut_pre$RTstatus)
       focal_mut_pre_mean_cp$unique_sample_id<-unique(focal_mut_pre$unique_sample_id.y)
       
       ### calulate mean of CP for postRT samples
       focal_mut_post<-mutaion_RT_focal[grepl("post",mutaion_RT_focal$identifier),]
       focal_mut_post_mean_cp<-aggregate( cellular_prevalence ~ new_cluster_id, focal_mut_post, mean )
       focal_mut_post_mean_cp$RTstatus<-unique(focal_mut_post$RTstatus)
       focal_mut_post_mean_cp$unique_sample_id<-unique(focal_mut_post$unique_sample_id.y)
      
       both<-rbind(focal_mut_pre_mean_cp,focal_mut_post_mean_cp)
       both$RTstatus<-gsub("no","before",both$RTstatus)
       out_res_py<-rbind(both, out_res_py)

}

#mutaion_RT<-merge(py_data_good,meta_paired, by.x = "sample_id", by.y = "sampleid")
#mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id.x,mutaion_RT$RTstatus, sep = "_")
#mutaion_RT$pos_identifier<-gsub("chr","",mutaion_RT$pos_id)
#samples<-unique(mutaion_RT$unique_sample_id.x)

paired_samples<-unique(out_res_py$unique_sample_id)
my.plot <- vector(mode = "list", length = length(paired_samples)) 


for (pp in 1:length(paired_samples)){

  both<-out_res_py[out_res_py$unique_sample_id == paired_samples[pp],]

  plot_line<-ggplot(both, aes(x=RTstatus, y=cellular_prevalence, group=as.factor(new_cluster_id))) +
  geom_line(aes(color=as.factor(new_cluster_id)))+
  ggtitle(paired_samples[pp]) +
  geom_point(aes(color=as.factor(new_cluster_id))) +
  xlab("Sample") + ylab("Cellular_Prevalence")+
  theme_bw() +
  guides(color=guide_legend("Cluster ID")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
        plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
        my.plot[[pp]]<-plot_line
}

pdf("~/Desktop/line_plot_paired_samples_good_clusters.pdf", width = 12, height = 22)
plots_line<-multiplot(plotlist = my.plot[1:7],cols= 2) 
print(plots_line)
dev.off()


####### scatter plot of CP ######
####### load mafs files #########
filenames_CN <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","AF")]
     })

shared_clonal_mutations<-NULL
for(ll in 1:length(attackStats_CN)){

    attackStats_CN_focal<-data.frame(attackStats_CN[ll])
    name<-unique(na.omit(attackStats_CN_focal)[,1])
    attackStats_CN_focal$sample_id[is.na(attackStats_CN_focal$sample_id)] <- name
    attackStats_CN_focal$Chromosome<-gsub("chr","",attackStats_CN_focal$Chromosome)
    shared_clonal_mutations<-rbind(attackStats_CN_focal,shared_clonal_mutations)
}



#mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange","RT_status","RT_staus_code")]
mutaion_maf<-merge(shared_clonal_mutations,meta_paired, by.x = "sample_id", by.y = "sampleid")
mutaion_maf<-mutaion_maf[mutaion_maf$Start!="777428",]
#unique(mutaion_RT$unique_sample_id.x)

out_res<-NULL
for(i in 1:length(paired_samples)){

   ## loop through each variant file
   paied_sample_focal<-mutaion_maf[mutaion_maf$unique_sample_id==paired_samples[i],]
   paied_sample_focal$mutation_id<-paste(paied_sample_focal$Chromosome,paied_sample_focal$Start , sep = "_")
   paied_sample_focal$vaf_corrected_purity<-(paied_sample_focal$AF *paied_sample_focal$purity)
   muttaions<-unique(paied_sample_focal$mutation_id)


      ## loop through pyclone file and adjust false cluster
       py_focal<-py_data_good[py_data_good$unique_sample_id==paired_samples[i],]
       ## remove a cluster with only 1 mutation and rename the clusters (new column "new_cluster_id" will be generated)##
       count_table<-table(py_focal$cluster_id)
       to_keep<-data.frame(which(count_table > length(unique(py_focal$sample_id))))
       py_focal_tokeep<-py_focal[py_focal$cluster_id%in%rownames(to_keep),]
       py_focal_tokeep$pos_id<-gsub("chr","",py_focal_tokeep$pos_id)

       py_focal_tokeep_renamed<-py_focal_tokeep |>
       arrange(cluster_id) |>
       mutate(new_cluster_id = cur_group_id()-1,
         .by = "cluster_id",
         .after = "cluster_id")

       for (j in 1:length(muttaions)){
        paied_sample_focal_mut<-paied_sample_focal[paied_sample_focal$mutation_id==muttaions[j],]
        paied_sample_focal_mut_pre<-paied_sample_focal_mut[grepl("no",paied_sample_focal_mut$RTstatus),]
        paied_sample_focal_mut_post<-paied_sample_focal_mut[grepl("post",paied_sample_focal_mut$RTstatus),]

        #mean_vaf_pre<-mean(paied_sample_focal_mut_pre$vaf_corrected_purity)
        mean_vaf_pre<-mean(paied_sample_focal_mut_pre$AF)

        vaf_table<-data.frame("vaf_pre"= mean_vaf_pre)
        #vaf_table$vaf_post<-mean(paied_sample_focal_mut_post$vaf_corrected_purity)
        vaf_table$vaf_post<-mean(paied_sample_focal_mut_post$AF)
        vaf_table$Gene<-unique(paied_sample_focal_mut$Gene)
        vaf_table$sample_status<-paired_samples[i]
        vaf_table$mutation<-muttaions[j]

        py_focal_tokeep_foc<-py_focal_tokeep_renamed[py_focal_tokeep_renamed$new_cluster_id%in%rownames(to_keep),]
        cluster<-py_focal_tokeep_foc[py_focal_tokeep_foc$pos_id==muttaions[j],]

        #pyclone_focal<-mutaion_RT[mutaion_RT$unique_sample_id.y== samples[i],]
        #cluster<-pyclone_focal[pyclone_focal$pos_identifier==muttaions[j],]
        if (length(cluster$new_cluster_id)>=1){         
          vaf_table$new_cluster_id<-unique(cluster$new_cluster_id)
      
        } 
        
        if (length(cluster$new_cluster_id)==0){
          vaf_table$new_cluster_id<-"not_available"
        }
        
        vaf_table
#       ccf_table$cluetr_id_pre<-unique(focal_mut_pre$cluster_id)
#       ccf_table$cluetr_id_post<-unique(focal_mut_post$cluster_id)
#
        out_res<-rbind(vaf_table,out_res)

  }
}

out_res[is.na(out_res)] <- 0

out_res$driver<-ifelse(out_res$Gene=="TP53","TP53",
                ifelse(out_res$Gene=="ATRX","ATRX",
                ifelse(out_res$Gene=="MUC16","MUC16",
                ifelse(out_res$Gene=="TNC","TNC",
                ifelse(out_res$Gene=="CLTC","CLTC",
                ifelse(out_res$Gene=="NBEA","NBEA",
                ifelse(out_res$Gene=="RBM10","RBM10",
                ifelse(out_res$Gene=="COL2A1","COL2A1",

                "-"))))))))

out_res
out_res<-out_res[out_res$new_cluster_id!="not_available",]

#out_res<-out_res[!grepl("not_available",out_res$cluster_id),]

#pdf(file = "~/Desktop/scatter_plot_VAF_corrected_purity_SRC167_merged_regions_.pdf", width = 7, height = )
#plotC<-ggplot(out_res, aes(x = vaf_pre, y = vaf_post, color = as.factor(cluster_id))) +
#geom_point(aes(color = as.factor(cluster_id)), size = 4,alpha = 0.6)
#geom_point(alpha= 0.7, size = 3) +
#geom_rug(position = "jitter", size = 0.7, color = "black")+
#xlab("VAF_preRT(corrected for purity)") + ylab("VAF_postRT(corrected for purity)")+
#  theme_bw() +
#  guides(color=guide_legend("Cluster ID")) +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"),
#        legend.key = element_rect(fill = NA, colour = NA, size = 0.25),axis.text=element_text(size=15),axis.title=element_text(size=15),
#        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))) +
#        geom_text_repel(data=subset(out_res, driver == "TP53" | driver == "ATRX"),
#        aes(vaf_pre,vaf_post,label=driver),color="black",size = 3.5)+
#         theme(legend.position = "bottom")
#print(plotC)
#dev.off()

### make a scatter plot and annotate with 
# "SRC125" "SRC127" "SRC167" "SRC169" "SRC170" "SRC171" "TB9051"

paired_samples_scatter<-unique(out_res$sample_status)
my.plot.scatter <- vector(mode = "list", length = length(paired_samples_scatter)) 

for (ss in 1:length(paired_samples_scatter)){

     out_res_good<-out_res[out_res$sample_status == paired_samples_scatter[ss],]
     plotC<-ggplot(out_res_good, aes(x = vaf_pre, y = vaf_post, color = as.factor(new_cluster_id))) +
     #geom_point(aes(color = as.factor(cluster_id)), size = 4,alpha = 0.6)
     ggtitle(paired_samples_scatter[ss]) +
     geom_point(alpha= 0.8, size = 3.5) +
     geom_rug(position = "jitter", size = 0.7, color = "black")+
     xlab("VAF_preRT(corrected for purity)") + ylab("VAF_postRT(corrected for purity)")+
     theme_bw() +
     guides(color=guide_legend("Cluster ID")) +
     geom_text_repel(data=subset(out_res_good, driver == "TP53" | driver == "ATRX"| driver == "TNC" | driver == "MUC16"| driver == "CLTC" | driver == "NBEA"| driver == "RBM10" | driver == "COL2A1"),
                     aes(vaf_pre,vaf_post,label=driver),color="black",size = 3.5)+
                     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = NA, colour = NA, size = 0.25),axis.text=element_text(size=15),axis.title=element_text(size=15),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))) 
     my.plot.scatter[[ss]]<-plotC
}

pdf("~/Desktop/scatter_plot_paired_samples_annotated_good_clusters.pdf", width = 15, height = 22)
plots_scatter<-multiplot(plotlist = my.plot.scatter[1:7],cols= 2) 
print(plots_scatter)
dev.off()


###########
### make heatmap of mutations across all regions
#https://rpubs.com/crazyhottommy/heatmap_demystified
mutaion_maf_all<-merge(shared_clonal_mutations,meta, by.x = "sample_id", by.y = "sampleid")
samples_all<-unique(mutaion_maf_all$sample_id)


for(ss in 1:length(samples_all)){

    focal_muttaion_tab<-mutaion_maf_all[mutaion_maf_all$unique_sample_id=="SRC125",]
    focal_muttaion_tab<-focal_muttaion_tab[,c(1,8)]
    #all_muts<-split(focal_muttaion_tab, focal_muttaion_tab$sample_id)

    all_genes<-data.frame("Gene"=unique(focal_muttaion_tab$Gene))

    subject<-unique(focal_muttaion_tab$sample_id)
    my_array<-array (NA, c (nrow (all_genes),(length(subject))))
    rownames(my_array)<-all_genes$Gene
  
    for(mm in 1:length(subject)){

      focal_subject<-focal_muttaion_tab[focal_muttaion_tab$sample_id==subject[mm],]
      focal_subject_gene<-merge(focal_subject,all_genes, by="Gene", all = T)

      focal_subject_gene_match <- focal_subject_gene[match(all_genes$Gene, focal_subject_gene$Gene),]
      focal_subject_gene_match[is.na(focal_subject_gene_match)]<-0
      focal_subject_gene_match$sample_id<-gsub(subject[mm],1,focal_subject_gene_match$sample_id)
      names(focal_subject_gene_match)[2]<-subject[mm]
      my_array[,mm]<-focal_subject_gene_match[,2]
    }

      my_array_df<-data.frame(my_array)
      my_array_df_good<-data.frame(lapply(my_array_df,as.numeric))     
      colnames(my_array_df_good)<-subject
      my_array_df_good$gene<-rownames(my_array_df)

      mut.tidy<-my_array_df_good %>% tidyr::gather(SRC, my_array_df_good, 1:7)
      mut.tidy$my_array_df_good<- factor(mut.tidy$my_array_df_good)
      colnames(mut.tidy)<-c("gene","sample","mutated")
      gg<- ggplot(mut.tidy, aes(x=sample, y=gene, fill=mutated)) + geom_tile(color="white", size=0.1)

}



### ***** ####



###########
### line plot of vaf dor each muttaion
### load maf file
#pre_maf<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/SRC125_1.filtered.variants.oxomerge.final.txt", header = TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","AF")]
#pre_maf$pos_id<-paste(pre_maf$Chromosome,pre_maf$Start, sep = "_")
#
#post_maf<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/SRC125_7.filtered.variants.oxomerge.final.txt", header = TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","AF")]
#post_maf$pos_id<-paste(post_maf$Chromosome,post_maf$Start, sep = "_")
#
### cluster file for pre and post 
#py_pre<-py_data_good[(py_data_good$sample_id%in%pre_maf$sample_id),]
#py_post<-py_data_good[(py_data_good$sample_id%in%post_maf$sample_id),]
#
###
#cluster_pre<-merge(py_pre,pre_maf, by="pos_id", all.x=T)[,c("pos_id","Chrom","Pos","Ref.x","Alt.x","sample_id.x","cluster_id" ,"cellular_prevalence","cellular_prevalence_std","cluster_assignment_prob","Gene","Impact","AF")]
#cluster_pre$RT_code<-rep("First_preRT")
#cluster_pre$AF[is.na(cluster_pre$AF)] <- 0
#
#
#cluster_post<-merge(py_post,post_maf, by="pos_id", all.x=T)[,c("pos_id","Chrom","Pos","Ref.x","Alt.x","sample_id.x","cluster_id" ,"cellular_prevalence","cellular_prevalence_std","cluster_assignment_prob","Gene","Impact","AF")]
#cluster_post$RT_code<-rep("Last_postRT")
#cluster_post$AF[is.na(cluster_post$AF)] <- 0
#
#
#both_Ts<-rbind(cluster_pre,cluster_post)
#colnames(both_Ts)[4]<-"populations"
#both_Ts$cluster_id<-paste("cluster",both_Ts$cluster_id, sep = "_")
#
#
#pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sample_SRC167B.pdf", width = 9, height =9 )
#plotC<-ggplot (data =both_Ts , aes(x = RT_code, y = AF, group = interaction(cluster_id, Pos) ,color = cluster_id)) + 
#  geom_line(aes(color = cluster_id), size=0.7, alpha=0.7)+
#  labs(title = "SRC167_5-SRC167_1")+  
#  geom_text(data = both_Ts %>% filter(RT_code == "First_preRT") %>% filter(Impact == "nonsynonymous SNV"), 
#            aes(label = Gene) , 
#            hjust = 1.5, 
#            size = 1.5) +
#  xlab("Sample") + ylab("VAF") +
#  guides(colour = guide_legend(override.aes = list(alpha = 3)))+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"),
#        legend.key = element_rect(fill = NA, colour = NA, size = 2),axis.text=element_text(size=16),axis.title=element_text(size=18),
#        plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
#print(plotC)
#dev.off() 

 
#plot_vaf(both_Ts,VAF_plot)


####
