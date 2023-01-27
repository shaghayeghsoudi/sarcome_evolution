#### pyclone outputs ####
rm(list = ls())

library(tidyr)
library(dplyr)
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/PycloneResults/probes_add500_purity_full")

py_data<-read.table(file ="SRC125.pyclone.results.tsv", header = TRUE)
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))


py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")

### load maf file
pre_maf<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/SRC125_1.filtered.variants.oxomerge.final.txt", header = TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","AF")]
pre_maf$pos_id<-paste(pre_maf$Chromosome,pre_maf$Start, sep = "_")

post_maf<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/SRC125_7.filtered.variants.oxomerge.final.txt", header = TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","AF")]
post_maf$pos_id<-paste(post_maf$Chromosome,post_maf$Start, sep = "_")

### cluster file for pre and post 
py_pre<-py_data_good[(py_data_good$sample_id%in%pre_maf$sample_id),]
py_post<-py_data_good[(py_data_good$sample_id%in%post_maf$sample_id),]

###
cluster_pre<-merge(py_pre,pre_maf, by="pos_id", all.x=T)[,c("pos_id","Chrom","Pos","Ref.x","Alt.x","sample_id.x","cluster_id" ,"cellular_prevalence","cellular_prevalence_std","cluster_assignment_prob","Gene","Impact","AF")]
cluster_pre$RT_code<-rep("First_preRT")
cluster_pre$AF[is.na(cluster_pre$AF)] <- 0


cluster_post<-merge(py_post,post_maf, by="pos_id", all.x=T)[,c("pos_id","Chrom","Pos","Ref.x","Alt.x","sample_id.x","cluster_id" ,"cellular_prevalence","cellular_prevalence_std","cluster_assignment_prob","Gene","Impact","AF")]
cluster_post$RT_code<-rep("Last_postRT")
cluster_post$AF[is.na(cluster_post$AF)] <- 0


both_Ts<-rbind(cluster_pre,cluster_post)
#colnames(both_Ts)[4]<-"populations"
both_Ts$cluster_id<-paste("cluster",both_Ts$cluster_id, sep = "_")


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sample_SRC167B.pdf", width = 9, height =9 )
plotC<-ggplot (data =both_Ts , aes(x = RT_code, y = AF, group = interaction(cluster_id, Pos) ,color = cluster_id)) + 
  geom_line(aes(color = cluster_id), size=0.7, alpha=0.7)+
  labs(title = "SRC167_5-SRC167_1")+  
  geom_text(data = both_Ts %>% filter(RT_code == "First_preRT") %>% filter(Impact == "nonsynonymous SNV"), 
            aes(label = Gene) , 
            hjust = 1.5, 
            size = 1.5) +
  xlab("Sample") + ylab("VAF") +
  guides(colour = guide_legend(override.aes = list(alpha = 3)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = NA, colour = NA, size = 2),axis.text=element_text(size=16),axis.title=element_text(size=18),
        plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
print(plotC)
dev.off() 

 
plot_vaf(both_Ts,VAF_plot)


####
