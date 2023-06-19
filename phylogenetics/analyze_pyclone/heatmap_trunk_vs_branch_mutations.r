### make heatmap of branch vs. trunk mutations and annotate driver genes
rm(list = ls())

### analyze trees

#### pyclone outputs ####
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)

###################################
### load tree Json files ###
###################################
## read tree location file
tree<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/PycloneResults/pyclone_multiple_Titan_deepsequenced_newsolutions/paired_cluster_location.txt", header = TRUE)
tree<-tree[tree$sample_id!="TB12052",]
samples<-unique(tree$sample_id)

###################################
### load pyclone cluster files ###
###################################
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/PycloneResults/pyclone_multiple_Titan_deepsequenced_newsolutions")
filenames_py <- list.files("00-inputs-pyclone-full",pattern="*.pyclone.results.tsv", full.names = TRUE)
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

py_data <- do.call("rbind", attackStats_py) 
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))

py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")
py_data_good$pos_id<-gsub("chr","",py_data_good$pos_id)

###################################
### load metadata ###
################################### 
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")

meta$RTstatus<-ifelse(meta$sequenceofsamplevRT== "beforeRT", "noRT",
                        ifelse(meta$sequenceofsamplevRT== "nopreopRT", "noRT",
                        ifelse(meta$sequenceofsamplevRT== "afterRT", "postRT",
                        "-")))

meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)
mutaion_RT<-merge(py_data_good,meta, by.x = "sample_id", by.y = "sampleid")

mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id.x,mutaion_RT$RTstatus, sep = "_")
mutaion_RT$pos_identifier<-gsub("chr","",mutaion_RT$pos_id)

###################################
### load maf/variant files ###
###################################
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
mutaion_maf<-merge(shared_clonal_mutations,meta, by.x = "sample_id", by.y = "sampleid")
mutaion_maf<-mutaion_maf[mutaion_maf$Start!="777428",]

### highly muttated genes
driver<-c("TP53","MUC16","TNC","CLTC","NBEA","ATRX","RBM10","COL2A1")

####################################################################
### make heatmap based on muttaions on trunk, pre, post or both ###
####################################################################

#out_res<-NULL
#samples
# [1] "SRC125"  "SRC127"  "SRC167"  "SRC169"  "SRC170"  "SRC171"  "TB9051" 
# [8] "SRC130"  "SRC150"  "SRC168"  "SRC173"  "TB8016"  "TB9573"  "TB11985"
#[15] "TB13092" "TB13712" "TB13959" "TB22446"

for(i in 1:length(samples)){  ### loop through each sample
    
    ## TO DO: 
    ### read relevant json file to extact trunk and branch information
    tree_focal_raw<-tree[tree$sample_id==samples[i],]
    
    ## focal pyclone
    pyclona_focal<-py_data_good[py_data_good$unique_sample_id==samples[i],]
    pyclona_focal<-merge(pyclona_focal,tree_focal_raw, by.x = "cluster_id", by.y = "cluster_id")

    uniq<-pyclona_focal[!duplicated(pyclona_focal$pos_id),]
    count_table<-table(uniq$cluster_id)
    to_keep<-data.frame(which(count_table > 1))
    pyclona_focal_tokeep<-pyclona_focal[pyclona_focal$cluster_id%in%rownames(to_keep),]


     pyclona_focal_tokeep_renamed<-pyclona_focal_tokeep |>
     arrange(cluster_id) |>
     mutate(new_cluster_id = cur_group_id()-1,
         .by = "cluster_id",
         .after = "cluster_id")

        
    pyclona_focal_tokeep_renamed$tree_location<-gsub('[0-9]+', '', pyclona_focal_tokeep_renamed$tree_location)
    pyclona_focal_tokeep_renamed$tree_location<-paste(pyclona_focal_tokeep_renamed$tree_location, pyclona_focal_tokeep_renamed$new_cluster_id, sep = "_")
    
    py_clusters<-unique(pyclona_focal_tokeep_renamed$tree_location_detail)

    ## load maf file for the looped sample
    paied_sample_focal<-mutaion_maf[mutaion_maf$unique_sample_id==samples[i],]
    paied_sample_focal$mutation_id<-paste(paied_sample_focal$Chromosome,paied_sample_focal$Start , sep = "_")
    paied_sample_focal<-paied_sample_focal[,c(1,8,18,19,20)]
    paied_sample_focal_uniq<-paied_sample_focal[!duplicated(paied_sample_focal$mutation_id),c(2:5)]
    
    all_genes<-data.frame("Gene"=unique(paied_sample_focal_uniq$Gene))  ### all genes present in the focal sample
    my_array<-array (NA, c (nrow (all_genes),length(py_clusters)))
    rownames(my_array)<-all_genes$Gene


    out_res_cluster<-NULL
    for(cc in 1:length(py_clusters)){   ### loop through trunk or each branch mutations

           #for(cc in 1:6){
              pyclona_cluster<-pyclona_focal_tokeep_renamed[pyclona_focal_tokeep_renamed$tree_location_detail==py_clusters[cc],]  ### 
              pyclona_cluster_uniq<-pyclona_cluster[!duplicated(pyclona_cluster$pos_id),]
              
              if (length(unique(pyclona_cluster_uniq$pos_id))>1){
              pyclone_gene<-merge(pyclona_cluster_uniq,paied_sample_focal_uniq, by.x = "pos_id",by.y = "mutation_id")[,c("pos_id","new_cluster_id","tree_location","unique_sample_id.x","tree_location_detail","Gene", "RTstatus")]
              #pyclone_gene$TB_info<-rep(tree_focal[tree_focal$clustre_id%in%pyclone_gene$cluster_id,2])
              pyclone_gene_good<-pyclone_gene[,c("tree_location_detail","Gene")]

              focal_subject_gene<-merge(pyclone_gene_good,all_genes, by="Gene", all = T)  ## merge with all genes in the sample
              focal_subject_gene_match <- focal_subject_gene[match(all_genes$Gene, focal_subject_gene$Gene),]
              focal_subject_gene_match$tree_location_detail[!is.na(focal_subject_gene_match$tree_location_detail)] <- 1
              focal_subject_gene_match[is.na(focal_subject_gene_match)]<-0
              names(focal_subject_gene_match)[2]<-py_clusters[cc]
              my_array[,cc]<-focal_subject_gene_match[,2]
              }
    }

            my_array_df<-data.frame(my_array)
            colnames(my_array_df)<-py_clusters   
            my_array_df<-my_array_df[,which(unlist(lapply(my_array_df, function(x) !all(is.na(x)))))] ## remove a column with all NA
            my_array_df_good<-data.frame(lapply(my_array_df,as.numeric))     ### make them numeric
            rownames(my_array_df_good)<-rownames(my_array_df)
            
            
            my_array_df_good_no0<-my_array_df_good[!apply(my_array_df_good, 1, function(x) all(x == 0)), ] ### remove genes not present on any tree position
            my_array_df_good_no0$gene<-rownames(my_array_df_good_no0)  ## genes whcih are present on trees (including trunk and branch)
            rownames(my_array_df_good_no0)<-NULL
            ####

             ### find order based on presence of genes on trunk, B1, B2, .. respectively
             out_res_order<-NULL
             for(uu in (length(my_array_df_good_no0)-1):1){

               my_orde<-my_array_df_good_no0[,c(uu,ncol(my_array_df_good_no0))]
               my_orde<-my_orde[grep (1,my_orde[,1]),]
               names(my_orde)<-c("present","gene")
               out_res_order<-rbind(my_orde,out_res_order)
            }
            
               order_fill<-unique(out_res_order$gene)
            
            # order based on presentce at trunk
            #aa<-my_array_df_good_no0[,c("trunk","gene")]
            #order_fill<-aa[order(aa$trunk, decreasing = TRUE),2]   ## ordered 
            my_array_df_good_no0_sort <- my_array_df_good_no0[match(order_fill, my_array_df_good_no0$gene),] ## ordered based on the presence of genes on trunk
      

           mut.tidy<-my_array_df_good_no0_sort %>% tidyr::gather(tree_location, my_array_df_good_no0_sort, 1:ncol(my_array_df_good_no0_sort)-1)
           #mut.tidy<-my_array_df_good_no0 %>% tidyr::gather(tree_location, my_array_df_good_no0, 1:ncol(my_array_df_good_no0)-1)

           colnames(mut.tidy)<-c("gene","tree_location","mutated")
           mut.tidy$mutated<- factor(mut.tidy$mutated)
           
          mut.tidy$tree_location<-gsub("RT","",mut.tidy$tree_location)

          ### arrange ordering genes on heatmap
           axis_labels<-which(mut.tidy$gene%in%driver) # show only driver
           qq<-mut.tidy %>%
           mutate(gene= fct_relevel(gene, order_fill))
           #ggsave(file = paste("~/Desktop/heatmap_tree_location_",samples[i],".pdf"), height = 4, width = 2)
           ggplot(data = qq, aes(x=tree_location, y=gene, fill=mutated)) + geom_tile(color="white", size=0.1) +
           scale_y_discrete(breaks = mut.tidy$gene[axis_labels])+
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
           ggtitle(samples[i]) 
           ggsave(file = paste("~/Desktop/heatmap_tree_location_",samples[i],".pdf"), height = 4, width = 2)
}

### make heatmaps based on TB_id ####
###################
driver<-c("TP53","MUC16","TNC","CLTC","NBEA","ATRX","RBM10","COL2A1")

out_res<-NULL
for(i in 1:length(samples)){  ### loop through each sample
    
    ## TO DO: 
    ### read relevant json file to extact trunk and branch information
    ## tem solution
    tree_focal_raw<-tree[tree$sample_id==samples[i],]
    names(tree_focal_raw)[3]<-"TB_id"

    # drop row with NA and rename the clusters
    tree_focal<-tree_focal_raw %>% drop_na(TB_id)
    tree_focal<-tree_focal |>
     arrange(cluster_id) |>
     mutate(new_cluster_id = cur_group_id()-1,
         .by = "cluster_id",
         .after = "cluster_id")
 
    ## focal pyclone
    pyclona_focal<-py_data_good[py_data_good$unique_sample_id==samples[i],]
    ## remove a cluster with only 1 mutation and rename the clusters (new column "new_cluster_id" will be generated)##
    uniq<-pyclona_focal[!duplicated(pyclona_focal$pos_id),]
    count_table<-table(uniq$cluster_id)
    to_keep<-data.frame(which(count_table > 1))
    pyclona_focal_tokeep<-pyclona_focal[pyclona_focal$cluster_id%in%rownames(to_keep),]


     pyclona_focal_tokeep_renamed<-pyclona_focal_tokeep |>
     arrange(cluster_id) |>
     mutate(new_cluster_id = cur_group_id()-1,
         .by = "cluster_id",
         .after = "cluster_id")

    
    ### assgin location of clusters on the tree
    pyclona_focal<-merge(pyclona_focal_tokeep_renamed,tree_focal, by.x = "new_cluster_id", by.y = "new_cluster_id")
    py_clusters<-unique(tree_focal$TB_id)

    ## load maf file for the looped sample
    paied_sample_focal<-mutaion_maf[mutaion_maf$unique_sample_id==samples[i],]
    paied_sample_focal$mutation_id<-paste(paied_sample_focal$Chromosome,paied_sample_focal$Start , sep = "_")
    paied_sample_focal<-paied_sample_focal[,c(1,8,18,19,20)]
    paied_sample_focal_uniq<-paied_sample_focal[!duplicated(paied_sample_focal$mutation_id),c(2:5)]
    
    all_genes<-data.frame("Gene"=unique(paied_sample_focal_uniq$Gene))  ### all genes present in the focal sample
    my_array<-array (NA, c (nrow (all_genes),length(py_clusters)))
    rownames(my_array)<-all_genes$Gene


    out_res_cluster<-NULL
    for(cc in 1:length(py_clusters)){   ### loop through trunk or each branch mutations

           #for(cc in 1:6){
              pyclona_cluster<-pyclona_focal[pyclona_focal$TB_id==py_clusters[cc],]  ### 
              pyclona_cluster_uniq<-pyclona_cluster[!duplicated(pyclona_cluster$pos_id),]
              
              if (length(unique(pyclona_cluster_uniq$pos_id))>1){
              pyclone_gene<-merge(pyclona_cluster_uniq,paied_sample_focal_uniq, by.x = "pos_id",by.y = "mutation_id")[,c("pos_id","new_cluster_id","TB_id","unique_sample_id.x","tree_location_detail","Gene", "RTstatus")]
              #pyclone_gene$TB_info<-rep(tree_focal[tree_focal$clustre_id%in%pyclone_gene$cluster_id,2])
              pyclone_gene_good<-pyclone_gene[,c("tree_location_detail","Gene")]

              focal_subject_gene<-merge(pyclone_gene_good,all_genes, by="Gene", all = T)  ## merge with all genes in the sample
              focal_subject_gene_match <- focal_subject_gene[match(all_genes$Gene, focal_subject_gene$Gene),]
              focal_subject_gene_match$tree_location_detail[!is.na(focal_subject_gene_match$tree_location_detail)] <- 1
              focal_subject_gene_match[is.na(focal_subject_gene_match)]<-0
              names(focal_subject_gene_match)[2]<-py_clusters[cc]
              my_array[,cc]<-focal_subject_gene_match[,2]
              }
    }

            my_array_df<-data.frame(my_array)
            colnames(my_array_df)<-py_clusters   
            my_array_df<-my_array_df[,which(unlist(lapply(my_array_df, function(x) !all(is.na(x)))))] ## remove a column with all NA
            my_array_df_good<-data.frame(lapply(my_array_df,as.numeric))     ### make them numeric
            rownames(my_array_df_good)<-rownames(my_array_df)
            
            
            my_array_df_good_no0<-my_array_df_good[!apply(my_array_df_good, 1, function(x) all(x == 0)), ] ### remove genes not present on any tree position
            my_array_df_good_no0$gene<-rownames(my_array_df_good_no0)  ## genes whcih are present on trees (including trunk and branch)
            rownames(my_array_df_good_no0)<-NULL
            ####

             ### find order based on presence of genes on trunk, B1, B2, .. respectively
             out_res_order<-NULL
             for(uu in (length(my_array_df_good_no0)-1):1){

               my_orde<-my_array_df_good_no0[,c(uu,ncol(my_array_df_good_no0))]
               my_orde<-my_orde[grep (1,my_orde[,1]),]
               names(my_orde)<-c("present","gene")
               out_res_order<-rbind(my_orde,out_res_order)
            }
            
               order_fill<-unique(out_res_order$gene)
            
            # order based on presentce at trunk
            #aa<-my_array_df_good_no0[,c("trunk","gene")]
            #order_fill<-aa[order(aa$trunk, decreasing = TRUE),2]   ## ordered 
            my_array_df_good_no0_sort <- my_array_df_good_no0[match(order_fill, my_array_df_good_no0$gene),] ## ordered based on the presence of genes on trunk
      

           mut.tidy<-my_array_df_good_no0_sort %>% tidyr::gather(tree_location, my_array_df_good_no0_sort, 1:ncol(my_array_df_good_no0_sort)-1)
           #mut.tidy<-my_array_df_good_no0 %>% tidyr::gather(tree_location, my_array_df_good_no0, 1:ncol(my_array_df_good_no0)-1)

           colnames(mut.tidy)<-c("gene","tree_location","mutated")
           mut.tidy$mutated<- factor(mut.tidy$mutated)
           
          
          ### arrange ordering genes on heatmap
           axis_labels<-which(mut.tidy$gene%in%driver) # show only driver
           qq<-mut.tidy %>%
           mutate(gene= fct_relevel(gene, order_fill))
           ggsave("xx.pdf", height = 4, width = 2)
           ggplot(data = qq, aes(x=tree_location, y=gene, fill=mutated)) + geom_tile(color="white", size=0.1) +
           scale_y_discrete(breaks = mut.tidy$gene[axis_labels])
}
