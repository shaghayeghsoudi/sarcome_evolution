#### pyclone outputs ####
rm(list = ls())

library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

###
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/PycloneResults/probes_add500_purity_full")
filenames_py <- list.files("paired_samples",pattern="*.pyclone.results.tsv", full.names = TRUE)
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

py_data <- do.call("rbind", attackStats_py) 
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))

py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")


## load meta data 
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)
paied_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167")
meta_paied<-meta[meta$unique_sample_id%in%paied_samples,]

py_samples<-unique(py_data_good$unique_sample_id)

out_res_py<-NULL
for (zz in 1:length(py_samples)){

       py_focal<-py_data_good[py_data_good$unique_sample_id==py_samples[zz],]
       mutaion_RT_focal<-merge(py_focal,meta_paied, by.x = "sample_id", by.y = "sample_id")
       mutaion_RT_focal$identifier<-paste(mutaion_RT_focal$unique_sample_id.x,mutaion_RT_focal$RT_status, sep = "_")
       mutaion_RT_focal$pos_identifier<-gsub("chr","",mutaion_RT_focal$pos_id)

       ### calulate mean of CP for preRT samples
       focal_mut_pre<-mutaion_RT_focal[grepl("pre",mutaion_RT_focal$identifier),]
       focal_mut_pre_mean_cp<-aggregate( cellular_prevalence ~ cluster_id, focal_mut_pre, mean )
       focal_mut_pre_mean_cp$RT_status<-unique(focal_mut_pre$RT_status)
       focal_mut_pre_mean_cp$unique_sample_id<-unique(focal_mut_pre$unique_sample_id.y)
       
       ### calulate mean of CP for postRT samples
       focal_mut_post<-mutaion_RT_focal[grepl("post",mutaion_RT_focal$identifier),]
       focal_mut_post_mean_cp<-aggregate( cellular_prevalence ~ cluster_id, focal_mut_post, mean )
       focal_mut_post_mean_cp$RT_status<-unique(focal_mut_post$RT_status)
       focal_mut_post_mean_cp$unique_sample_id<-unique(focal_mut_post$unique_sample_id.y)
      
       both<-rbind(focal_mut_pre_mean_cp,focal_mut_post_mean_cp)
       both$RT_status<-gsub("pre","before",both$RT_status)
       out_res_py<-rbind(both, out_res_py)

}



##################


#py_data<-read.table(file ="SRC167.pyclone.results.tsv", header = TRUE)
#py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))

#py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
#py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")

#meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
#colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
#meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)
#paied_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167")
#meta_paied<-meta[meta$unique_sample_id%in%paied_samples,]

#mutaion_RT<-merge(py_data_good,meta_paied, by.x = "sample_id", by.y = "sample_id")

#mutaion_RT$pos_id<-paste(mutaion_RT$Chrom,mutaion_RT$Pos, sep = "_")
#mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id.x,mutaion_RT$RT_status, sep = "_")
#mutaion_RT$pos_identifier<-gsub("chr","",mutaion_RT$pos_id)

#samples<-unique(mutaion_RT$unique_sample_id.x)

### scatter plot of CP ####
#out_res<-NULL
#for (i in 1:length(samples)){
#
#    mutaion_RT_focal<-mutaion_RT[mutaion_RT$unique_sample_id.x==samples[i],]
#    mutations<-unique(mutaion_RT_focal$pos_id)
#    
#       for (j in 1:length(mutations)){
#        #for (j in 1:50){
#
#
#       focal_mut<-mutaion_RT_focal[mutaion_RT_focal$pos_id==mutations[j],]
#       focal_mut_pre<-focal_mut[grepl("pre",focal_mut$identifier),]
#       
#       focal_mut_post<-focal_mut[grepl("post",focal_mut$identifier),]
#       
#       mean_ccf_pre<-mean(focal_mut_pre$cellular_prevalence)
#       ccf_table<-data.frame("ccf_pre"= mean_ccf_pre)
#       ccf_table$ccf_post<-mean(focal_mut_post$cellular_prevalence)
#
#       ccf_table$sample_status<-samples[i]
#       ccf_table$mutation<-mutations[j]
#
#       ccf_table$cluetr_id_pre<-unique(focal_mut_pre$cluster_id)
#       ccf_table$cluetr_id_post<-unique(focal_mut_post$cluster_id)
#
#       out_res<-rbind(ccf_table,out_res)
#    }  
#
##}

### scatter plot of Celleular prevalence ####
#ggplot(out_res, aes(x = ccf_pre, y = ccf_pre, color = as.factor(cluetr_id_pre))) +
#geom_point()




### line plot of CP per cluer ###

out_res<-NULL
for (i in 1:length(samples)){

       mutaion_RT_focal<-mutaion_RT[mutaion_RT$unique_sample_id.x==samples[i],]
    
       focal_mut_pre<-mutaion_RT_focal[grepl("pre",mutaion_RT_focal$identifier),]
       focal_mut_pre_mean_cp<-aggregate( cellular_prevalence ~ cluster_id, focal_mut_pre, mean )
       focal_mut_pre_mean_cp$RT_status<-unique(focal_mut_pre$RT_status)
       focal_mut_pre_mean_cp$unique_sample_id<-unique(focal_mut_pre$unique_sample_id.y)
       
       
       focal_mut_post<-mutaion_RT_focal[grepl("post",mutaion_RT_focal$identifier),]
       focal_mut_post_mean_cp<-aggregate( cellular_prevalence ~ cluster_id, focal_mut_post, mean )
       focal_mut_post_mean_cp$RT_status<-unique(focal_mut_post$RT_status)
       focal_mut_post_mean_cp$unique_sample_id<-unique(focal_mut_post$unique_sample_id.y)
      
       both<-rbind(focal_mut_pre_mean_cp,focal_mut_post_mean_cp)
       both$RT_status<-gsub("pre","before",both$RT_status)
       #out_res<-rbind(ccf_table,out_res)
 }  


p<-ggplot(both, aes(x=RT_status, y=cellular_prevalence, group=as.factor(cluster_id))) +
  geom_line(aes(color=as.factor(cluster_id)))+
  geom_point(aes(color=as.factor(cluster_id)))
p



pdf(file = "~/Desktop/line_plot_SRC167_merged_regions_CP.pdf", width = 8, height =8 )
plotB<-ggplot(data = both, aes(x=RT_status, y=cellular_prevalence, group=as.factor(cluster_id))) +
  geom_line(aes(color = as.factor(cluster_id)), size = 2,alpha = 0.6) +
  #labs(title = paste("Best Tree",gsub(".json", "",best_tree_fileID), sep = " "))+  
  geom_point(aes(color = as.factor(cluster_id)), size = 4,alpha = 0.6) +
  #  Labelling as desired
  xlab("Sample") + ylab("Cellular_Prevalence")+
  theme_bw() +
  guides(color=guide_legend("Cluster ID")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = NA, colour = NA, size = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
        plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
        print(plotB)
        dev.off() 


    
###################################
###################################
### scatter plot of VAF corrected for purity

#setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/")

## load mafs
filenames_CN <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.csv(x, header=TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","AF")]
     })


shared_clonal_mutations <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix<-shared_clonal_mutations%>%drop_na(sample_id)
shared_clonal_mutations_fix$Chromosome<-gsub("chr","",shared_clonal_mutations_fix$Chromosome)

#mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange","RT_status","RT_staus_code")]

mutaion_maf<-merge(shared_clonal_mutations_fix,meta_paied, by.x = "sample_id", by.y = "sample_id")
mutaion_maf<-mutaion_maf[mutaion_maf$Start!="777428",]


mutaion_maf$mutaion_maf<-gsub("noRT-1","noRT",mutaion_maf$RT_status)
mutaion_maf$mutaion_maf<-gsub("noRT-2","noRT",mutaion_maf$RT_status)

mutaion_maf$RTtretment<-ifelse(mutaion_maf$RT_status== "preRT", "naive",
                        ifelse(mutaion_maf$RT_status== "noRT", "naive",
                        ifelse(mutaion_maf$RT_status== "postRT", "treatment",
                        "-")))


#unique(mutaion_RT$unique_sample_id.x)
samples<-"SRC167"

out_res<-NULL
for(i in 1:length(samples)){

   paied_sample_focal<-mutaion_maf[mutaion_maf$unique_sample_id==samples[i],]
   paied_sample_focal$mutation_id<-paste(paied_sample_focal$Chromosome,paied_sample_focal$Start , sep = "_")
   paied_sample_focal$vaf_corrected_purity<-(paied_sample_focal$AF *paied_sample_focal$purity)
   muttaions<-unique(paied_sample_focal$mutation_id)

  for (j in 1:length(muttaions)){
    #for (j in 1:20){

        paied_sample_focal_mut<-paied_sample_focal[paied_sample_focal$mutation_id==muttaions[j],]
        paied_sample_focal_mut_pre<-paied_sample_focal_mut[grepl("pre",paied_sample_focal_mut$RT_status),]
        paied_sample_focal_mut_post<-paied_sample_focal_mut[grepl("post",paied_sample_focal_mut$RT_status),]

        mean_vaf_pre<-mean(paied_sample_focal_mut_pre$vaf_corrected_purity)
        vaf_table<-data.frame("vaf_pre"= mean_vaf_pre)
        vaf_table$vaf_post<-mean(paied_sample_focal_mut_post$vaf_corrected_purity)
        
        vaf_table$sample_status<-samples[i]
        vaf_table$mutation<-muttaions[j]

        cluster<-mutaion_RT[mutaion_RT$pos_identifier==muttaions[j],]
        if (length(cluster$cluster_id)>=1){
          
          vaf_table$cluster_id<-unique(cluster$cluster_id)
        

        } else {
          vaf_table$cluster_id<-"not_available"
        }
        
#       ccf_table$cluetr_id_pre<-unique(focal_mut_pre$cluster_id)
#       ccf_table$cluetr_id_post<-unique(focal_mut_post$cluster_id)
#
        out_res<-rbind(vaf_table,out_res)


  }
}


out_res[is.na(out_res)] <- 0
out_res<-out_res[!grepl("not_available",out_res$cluster_id),]

pdf(file = "~/Desktop/scatter_plot_VAF_corrected_purity_SRC167_merged_regions_.pdf", width = 7, height = )
plotC<-ggplot(out_res, aes(x = vaf_pre, y = vaf_post, color = as.factor(cluster_id))) +
#geom_point(aes(color = as.factor(cluster_id)), size = 4,alpha = 0.6)
geom_point(alpha= 0.7, size = 3) +
geom_rug(position = "jitter", size = 0.5, color = "black")+
xlab("VAF_preRT(corrected for purity)") + ylab("VAF_postRT(corrected for purity)")+
  theme_bw() +
  guides(color=guide_legend("Cluster ID")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = NA, colour = NA, size = 0.25),axis.text=element_text(size=15),axis.title=element_text(size=15),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))

print(plotC)
dev.off()


##### end #######

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
