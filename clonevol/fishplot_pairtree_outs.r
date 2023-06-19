


#### make fish tree plots for pyclone outputs
rm(list = ls())
#library(fishplot)
library(tidyverse)
library(clonevol)
#library("plyr")
library("rjson")

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/trees_from_updated_solutions/data_and_trees/")
####################
## load meta data ##
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")
#colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)
#paied_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167")
#meta_paied<-meta[meta$unique_sample_id%in%paied_samples,]
#py_samples<-unique(py_data_good$unique_sample_id)


filenames_py <- list.files("00-inputs-pyclone",pattern="*.pyclone.results.tsv", full.names = TRUE)  ### load pyclone cluster files
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

py_data <- do.call("rbind", attackStats_py) 
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))
#py_data_good<-py_data_good[!grepl("_C",py_data_good$sample_id),]
py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")

### load variant files 
vafs<-list.files ("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants", pattern = "*.filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_mafs <- lapply(vafs,function(x) {
     read.delim(x, header=TRUE, sep = "\t")[,c(1:7)]
     })

for (i in 1:length(attackStats_mafs)){
    attackStats_mafs[[i]]<-cbind(attackStats_mafs[[i]],vafs[i])
    }
variants_aa <- do.call("rbind", attackStats_mafs) 
variants_aa$sample_id<-gsub("/Users/shsoudi/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/", "",gsub(".filtered.variants.oxomerge.final.txt","",variants_aa[,8]))
variants_aa$uniq_sample_id<-gsub("_.*$","",variants_aa$sample_id)

drivers<-c("TP53","ATRX","MUC16","TNC","CLTC","NBEA","RBM10","COL2A1")

### create additional columns to identify drivers
variants_aa$gene<-ifelse(variants_aa$Gene== "TP53", "TP53",              
                        ifelse(variants_aa$Gene== "ATRX", "ATRX",
                        ifelse(variants_aa$Gene== "MUC16", "MUC16",
                        ifelse(variants_aa$Gene== "TNC", "TNC",
                        ifelse(variants_aa$Gene== "CLTC", "CLTC",
                        ifelse(variants_aa$Gene== "NBEA", "NBEA",
                        ifelse(variants_aa$Gene== "RBM10", "RBM10",
                        ifelse(variants_aa$Gene== "COL2A1", "COL2A1",
                        "-"))))))))

variants_aa$is.driver<-ifelse(variants_aa$Gene== "TP53", "TP53",
                        ifelse(variants_aa$Gene== "ATRX", "ATRX",
                        ifelse(variants_aa$Gene== "MUC16", "MUC16",
                        ifelse(variants_aa$Gene== "TNC", "TNC",
                        ifelse(variants_aa$Gene== "CLTC", "CLTC",
                        ifelse(variants_aa$Gene== "NBEA", "NBEA",
                        ifelse(variants_aa$Gene== "RBM10", "RBM10",
                        ifelse(variants_aa$Gene== "COL2A1", "COL2A1",
                        "FALSE"))))))))

variants_aa$pos_id<-paste(variants_aa$Chromosome ,variants_aa$Start ,sep = "_")
variants_aa$pos_id<-gsub("chr","",variants_aa$pos_id)

subject_id<-unique(py_data_good$unique_sample_id)


## create a new cluetr 0 , pyclone file format (normal cell with all cff= 1 based on pairtree)
out_res_var<-NULL
for(ss in 1:length(subject_id)){
  #out_res_all<-NULL
    
    ## remove cluster with only one mutation and rename the clusters based on the removed cluster 
    py_focal<-py_data_good[py_data_good$unique_sample_id==subject_id[ss],]
    uniq<-py_focal[!duplicated(py_focal$pos_id),]
    count_table<-table(uniq$cluster_id)
    to_keep<-data.frame(which(count_table > 1))
    py_focal_tokeep<-py_focal[py_focal$cluster_id%in%rownames(to_keep),]   ### keep clusters with more than one mutation


    py_focal_tokeep_renamed<-py_focal_tokeep |>    ### rename clusters based on the removed cluster
     arrange(cluster_id) |>
     mutate(new_cluster_id = cur_group_id()-1,
         .by = "cluster_id",
         .after = "cluster_id")
    
    py_focal_tokeep_renamed<-py_focal_tokeep_renamed[ , -which(names(py_focal_tokeep_renamed) %in% c("cluster_id"))]
    
    py_focal_tokeep_renamed$cluster_id<-(py_focal_tokeep_renamed$new_cluster_id)+1
    py_focal_tokeep_renamed$new_cluster_id<-py_focal_tokeep_renamed$cluster_id  ## replace new cluster ID column with cluster ID
    py_focal_tokeep_renamed<-py_focal_tokeep_renamed[1:(length(py_focal_tokeep_renamed)-1)]
    names(py_focal_tokeep_renamed)[6]<-"cluster_id"


    var_focal<-aa[aa$uniq_sample_id==subject_id[ss],]   ### variant file 
    not_inpy<-var_focal[!(var_focal$pos_id%in%py_focal$pos_id),]

    not_inpy$cluster_id<-0
    not_inpy$cellular_prevalence<-1
    not_inpy$cellular_prevalence_std<-0
    not_inpy$cluster_assignment_prob<-1
    cluster_0<-unique(not_inpy$pos_id)

    out_res_pos<-NULL
    for(rr in 1:length(cluster_0)){  

        sudo_positions_focal<-not_inpy[not_inpy$pos_id==cluster_0[rr],]
        sudo_positions_focal<-sudo_positions_focal[1,]
        sudo_positions_focal<-sudo_positions_focal[rep(seq_len(nrow(sudo_positions_focal)), each = length(unique(py_focal_tokeep_renamed$sample_id))), ]   ## repeat each row n= number of samples time 
        sudo_positions_focal<-sudo_positions_focal[,c("Chromosome","Start","Ref","Alt","uniq_sample_id","pos_id", "cluster_id", "cellular_prevalence", "cellular_prevalence_std","cluster_assignment_prob")]
        a<-rep(1:length(unique(py_focal_tokeep_renamed$sample_id)))
        sudo_positions_focal<-cbind(sudo_positions_focal,a)
        sudo_positions_focal$a<- paste(subject_id[ss], sudo_positions_focal$a, sep="_")
        sudo_positions_focal<-sudo_positions_focal[,c("Chromosome","Start","Ref","Alt","a","cluster_id","cellular_prevalence","cellular_prevalence_std","cluster_assignment_prob","uniq_sample_id","pos_id")]

        colnames(sudo_positions_focal)[c(1,2,5,10)]<-c("Chrom","Pos","sample_id","unique_sample_id")
        out_res_pos<-rbind(sudo_positions_focal,out_res_pos)        
    }

     out_res_var<-rbind(out_res_pos,py_focal_tokeep_renamed)
     out_res_var$cluster_id<-(out_res_var$cluster_id)+1
}


################
out_res_clon<-NULL
for (zz in 1:length(subject_id)){

    py_focal<-out_res_var[out_res_var$unique_sample_id==subject_id[zz],]
    mutaion_RT_focal<-merge(py_focal,meta, by.x = "sample_id", by.y = "sampleid")

    mutaion_RT_focal$RT_status<-ifelse(mutaion_RT_focal$sequenceofsamplevRT== "beforeRT", "preRT",
                        ifelse(mutaion_RT_focal$sequenceofsamplevRT== "nopreopRT", "preRT",
                        ifelse(mutaion_RT_focal$sequenceofsamplevRT== "afterRT", "postRT",
                        "-")))
    mutaion_RT_focal$identifier<-paste(mutaion_RT_focal$unique_sample_id.x,mutaion_RT_focal$RT_status, sep = "_")
    mutaion_RT_focal$code<-paste(mutaion_RT_focal$sample_id,mutaion_RT_focal$RT_status, sep = "_")
    mutaion_RT_focal$pos_identifier<-gsub("chr","",mutaion_RT_focal$pos_id)

    pos_identifier<-unique(mutaion_RT_focal$pos_identifier)
    mutaion_RT_focal_good<-mutaion_RT_focal[,c("code","cluster_id","cellular_prevalence","pos_identifier","sample_id")]
    
     out_res_phy<-NULL
     for(ii in 1:length(pos_identifier)) {

    focal_pos<-mutaion_RT_focal_good[mutaion_RT_focal_good$pos_identifier== pos_identifier[ii],]
    focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)*100
    #focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)/2


    #focal_pos<-out_res_mutaion_RT_focal_good[out_res_mutaion_RT_focal_good$pos_identifier== pos_identifier[ii],]
    #focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)*100

    aa<-data.frame(t(focal_pos))
    colnames(aa)<-aa[1,]
    bb<-aa[3,]
    rownames(bb)<-NULL
    colnames(bb)<-paste(colnames(bb),"ccf", sep = ".")
    bb[] <- lapply(bb, function(x) as.numeric(as.character(x)))
        
    vaf<-bb/2
    colnames(vaf)<-gsub("ccf","vaf",colnames(vaf))
    
    vaf2<-vaf
    colnames(vaf2)<-gsub(".vaf","",colnames(vaf2))
    
    both<-cbind(vaf,bb,vaf2)
    both$cluster<-unique(focal_pos$cluster_id)
    both$position<-unique(focal_pos$pos_identifier)
    variant_focal<-variants_aa[variants_aa$uniq_sample_id ==subject_id[zz] & variants_aa$pos_id==pos_identifier[ii],c("gene","is.driver")]
    both_variant_focal<-cbind(variant_focal,both)
  
    out_res_phy<-rbind(bb,out_res_phy)
    write.table(file = "out_res_phy", file = paste("~/Desktop/out_res_variants_clonevol_format_",subject_id[zz],".txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
       
       } ## pos_identifier loop
} ##

###################################################
#### take json files to prepare tree.tsv files ####
###################################################

file_trees<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/pairtree/paired_pyclone/jsons", pattern = "*.tree.json", full.names = TRUE)
json_file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/pairtree/paired_pyclone/jsons/SRC125.tree.json"
json_data <- fromJSON(paste(readLines(json_file), collapse=""))

 for (ii in 1:length(file_trees)){
     
     focal_file_tree<-file_trees[ii] 
     json_focal <- fromJSON(paste(readLines(focal_file_tree), collapse=""))  ## parse Json tree file into R

     ## adjust samples name to match with variant file
     subject<-gsub(".tree.json", "",sub(".*/", "", focal_file_tree)) 
     sample_id<-json_data$samples
     samples<-paste(subject,sample_id, sep = "_")
     samples_right<-gsub("R_","",samples)

     ### adjust clones and parent
     tree_focal<-data.frame("clone"=1:length(json_focal$phi))

     tree_focal<-data.frame("clone"=1:length(json_focal$phi))  ### howmany clones we get, including the normal cell clone
     focal_parent<-data.frame("parent"=json_focal$parents)
     focal_parent$parent<-(focal_parent$parent)+1
     data <- rbind(data.frame(parent = "-1"),focal_parent)
     data_clon_parent<-cbind(tree_focal,data)

     out_res_cell<-NULL
     for(kk in 1:length(json_data$phi)){

          focal<-json_data$phi[[kk]]*100
          decima<-as.numeric(format(round(focal,digits=2),nsmall=2))
          aa<-paste(as.character(decima), "%(0.95)/p=000", sep = "")
          bb<-paste(as.character(decima),aa, sep = "-")
          cc<-paste(samples_right,bb, sep = ":")
          cff_0<-cc[grep("0-0%",cc)]  ### remove clones with 0 CCF
          cc_non0<-cc[!(cc%in%cff_0)]
          dd<-data.frame("sample.with.nonzero.cell.frac.ci"=paste(cc_non0,collapse=","))
          dd$clone_id<-kk
          
          out_res_cell<-rbind(dd,out_res_cell)
          out_res_cell<-out_res_cell[sort(out_res_cell$clone_id, decreasing = TRUE),]

     }

}



