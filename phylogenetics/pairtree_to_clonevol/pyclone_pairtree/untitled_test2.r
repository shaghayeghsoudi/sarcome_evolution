
### full data set

#### make fish tree plots for pyclone outputs
rm(list = ls())
#library(fishplot)
library(tidyverse)
library(clonevol)
library("rjson")

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/trees_from_updated_solutions/data_and_trees/")

################################
######## load metadata #########
################################
meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated_shaghayegh_Feb2023_based_on_new_solutions.txt", header = TRUE, sep= "\t")
meta$unique_sample_id<-gsub("_.*$","",meta$sampleid)


meta$sequenceofsamplevRT <- gsub('beforeRT', 'preRT',
           gsub('afterRT', 'postRT',
           gsub('nopreopRT', 'noRT', meta$sequenceofsamplevRT )))


meta$sample_rt<-paste(meta$unique_sample_id,meta$sequenceofsamplevRT, sep = "_")
meta$sample_id<-paste(meta$sampleid , meta$sequenceofsamplevRT , sep = "_")
meta$subject_id<-paste(meta$unique_sample_id , meta$sequenceofsamplevRT , sep = "_")
#paired_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167","TB9052")
#paired_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167","TB9052")


###################################################
#### take json files to prepare tree.tsv files ####
###################################################

#file_trees<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/pairtree/paired_pyclone/jsons", pattern = "*.tree.json", full.names = TRUE)
#json_file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/pairtree/paired_pyclone/jsons/SRC125.tree.json"
#json_data <- fromJSON(paste(readLines(json_file), collapse=""))

### start from pairtree XX.tree.json file
file_trees<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/jsons", pattern = "*.tree.json", full.names = TRUE)
#json_data <- fromJSON(paste(readLines(file_trees), collapse=""))

## make subject IDs based on the samles you want to run
subject_pre<-sub('.*\\/', '', file_trees)
subject_id<-gsub(".tree.json","",subject_pre)
subject_id<-gsub("-","_",subject_id)   # "SRC125_postRT"


for (ii in 1:length(file_trees)){
     
     focal_file_tree<-file_trees[ii] 
     json_focal <- fromJSON(paste(readLines(focal_file_tree), collapse=""))  ## parse Json tree file into R
     ## adjust samples name to match with variant file
     subject_id<-gsub(".tree.json", "",sub(".*/", "", focal_file_tree)) 
     subject_id<-gsub("-","_",subject_id)
     
     subjecta<-gsub("_.*$","",subject_id)
     
     #subject<-gsub("-.*$","",subject)
     sample_id<-json_focal$samples
     sample_id<-gsub("R_","",sample_id)
     samples_right<-paste(subjecta,sample_id, sep = "_")

     ### adjust clones and parent
     tree_focal<-data.frame("clone"=1:length(json_focal$phi))  ### howmany clones we get, including the normal cell clone
     focal_parent<-data.frame("parent"=json_focal$parents)
     focal_parent$parent<-(focal_parent$parent)+1
     data <- rbind(data.frame(parent = "-1"),focal_parent)
     data_clon_parent<-cbind(tree_focal,data)

     out_res_cell<-NULL
     for(kk in 1:length(json_focal$phi)){

          focal<-json_focal$phi[[kk]]*100
          decima<-as.numeric(format(round(focal,digits=2),nsmall=2))
          aa<-paste(as.character(decima), "%(0.95)/p=000", sep = "")
          bb<-paste(as.character(decima),aa, sep = "-")
          cc<-paste(samples_right,bb, sep = ":")
          cff_0<-cc[grep("0-0%",cc)]  ### remove clones with 0 CCF
          cc_non0<-cc[!(cc%in%cff_0)]
          dd<-data.frame("sample.with.nonzero.cell.frac.ci"=paste(cc_non0,collapse=","))
          dd$clone_id<-kk
          
          dd_with_parent<-merge(dd,data_clon_parent, by.x = "clone_id", by.y= "clone")
          dd_with_parent<-dd_with_parent[,c("clone_id","parent","sample.with.nonzero.cell.frac.ci")]
          colnames(dd_with_parent)<-c("clone","parent","sample.with.nonzero.cell.frac.ci")

          out_res_cell<-rbind(dd_with_parent,out_res_cell)

     }

      write.table(out_res_cell, file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/tree.",subject_id,".tsv", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

}

################################
## load pyclone cluster files ##
################################
filenames_py <- list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/pyclone_newsolutions_for_pairtree_analysis/pyclone_tsv",pattern="*.pyclone.results.tsv", full.names = TRUE)  ### load pyclone cluster files
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

py_data <- do.call("rbind", attackStats_py) 

py_data_good <-merge(py_data,meta, by.x = "sample_id",by.y = "sampleid")[c("mutation_id","sample_id", "cluster_id", "cellular_prevalence","cellular_prevalence_std","cluster_assignment_prob","sample_rt","sample_id.y")]
py_data_good<-py_data_good %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))  ## make seperate columns for chrom, pos, ref and alt allels
#py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
#py_data_good$subject_id<-gsub("_.*$","",py_data_good$sample_id)

colnames(py_data_good)[colnames(py_data_good)=="sample_rt"]<-"subject_id"
colnames(py_data_good)[colnames(py_data_good)=="sample_id.y"]<-"sample_rt"
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")
py_data_good$pos_id<-gsub("chr","",py_data_good$pos_id)

py_data_good<-py_data_good[py_data_good$subject_id%in%subject_id,]


#py_data_good$RTstatus<-paste(py_data_good$sample_id,py_data_good,)


################################
##### load variant files #######
################################ 
vafs<-list.files ("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants", pattern = "*.filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_mafs <- lapply(vafs,function(x) {
     read.delim(x, header=TRUE, sep = "\t")[,c(1:7)]
     })

for (i in 1:length(attackStats_mafs)){
    attackStats_mafs[[i]]<-cbind(attackStats_mafs[[i]],vafs[i])
    }
variants <- do.call("rbind", attackStats_mafs) 
variants$sample_id<-gsub("/Users/shsoudi/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/filtered_variants/", "",gsub(".filtered.variants.oxomerge.final.txt","",variants[,8]))
#variants_aa$uniq_sample_id<-gsub("_.*$","",variants_aa$sample_id)

variants_aa<-merge(variants,meta, by.x = "sample_id", by.y ="sampleid" )[,c("Chromosome" , "Start" , "End" ,"Ref" ,"Alt" ,"Gene" ,"sample_id" ,"subject_id","sample_id.y")]
colnames(variants_aa)[colnames(variants_aa)=="sample_id.y"]<-"sample_rt"

#variants_aa$subject_id<-gsub("_.*$","",variants_aa$sample_id)
variants_aa<-variants_aa[variants_aa$subject_id%in%subject_id,]

drivers<-c("TP53","ATRX","MUC16","TNC","CLTC","NBEA","RBM10","COL2A1")

### create additional columns to specify when drivers are present, and if present make a sepereate column to mark them as "TRUE"
### This is required by clonevol
variants_aa$gene<-ifelse(variants_aa$Gene== "TP53", "TP53",              
                        ifelse(variants_aa$Gene== "ATRX", "ATRX",
                        ifelse(variants_aa$Gene== "MUC16", "MUC16",
                        ifelse(variants_aa$Gene== "TNC", "TNC",
                        ifelse(variants_aa$Gene== "CLTC", "CLTC",
                        ifelse(variants_aa$Gene== "NBEA", "NBEA",
                        ifelse(variants_aa$Gene== "RBM10", "RBM10",
                        ifelse(variants_aa$Gene== "COL2A1", "COL2A1",
                        "-"))))))))

variants_aa$is.driver<-ifelse(variants_aa$Gene== "TP53", "TRUE",
                        ifelse(variants_aa$Gene== "ATRX", "TRUE",
                        ifelse(variants_aa$Gene== "MUC16", "TRUE",
                        ifelse(variants_aa$Gene== "TNC", "TRUE",
                        ifelse(variants_aa$Gene== "CLTC", "TRUE",
                        ifelse(variants_aa$Gene== "NBEA", "TRUE",
                        ifelse(variants_aa$Gene== "RBM10", "TRUE",
                        ifelse(variants_aa$Gene== "COL2A1", "TRUE",
                        "FALSE"))))))))

variants_aa$pos_id<-paste(variants_aa$Chromosome ,variants_aa$Start ,sep = "_")
variants_aa$pos_id<-gsub("chr","",variants_aa$pos_id)

#subject_id<-unique(py_data_good$subject_id)
#subject_id<-paired_samples

#####################################################################################################################
#####################################################################################################################
## create a fake 0 and assign sudo variants to it based on pairtree for each sample and add to to cluster adjusted pyclone
#####################################################################################################################
#####################################################################################################################
#names(py_data_good)[5]<-"sample"
#names(py_data_good)[11]<-"sample_id"

out_res_var<-NULL
for(ss in 1:length(subject_id)){ 

#for(ss in 1:20){ 
   ### loop through each subject
  #out_res_all<-NULL
  #out_res_var<-NULL
    
    ## remove cluster with only one mutation and rename remaining clusters based on the removed cluster 
    py_focal<-py_data_good[py_data_good$subject_id==subject_id[ss],]
    #names(py_focal)[5]<-"sample"
    #names(py_focal)[11]<-"sample_id"

    
    
    uniq<-py_focal[!duplicated(py_focal$pos_id),]
    count_table<-table(uniq$cluster_id)
    to_keep<-data.frame(which(count_table > 1))
    py_focal_tokeep<-py_focal[py_focal$cluster_id%in%rownames(to_keep),]   ### keep clusters with more than one mutation


    py_focal_tokeep_renamed<-py_focal_tokeep |>    ### rename clusters based on the removed cluster, a "new_cluster_id" column will be generated 
     arrange(cluster_id) |>
     mutate(new_cluster_id = cur_group_id()-1,
         .by = "cluster_id",
         .after = "cluster_id")
    
    py_focal_tokeep_renamed<-py_focal_tokeep_renamed[ , -which(names(py_focal_tokeep_renamed) %in% c("cluster_id"))]
    
    py_focal_tokeep_renamed$cluster_id<-(py_focal_tokeep_renamed$new_cluster_id)+1   ### add 1 to each cluster to create fake cluster 0
    py_focal_tokeep_renamed$new_cluster_id<-py_focal_tokeep_renamed$cluster_id       ## replace new_cluster_ID column with cluster_ID column
    py_focal_tokeep_renamed<-py_focal_tokeep_renamed[1:(length(py_focal_tokeep_renamed)-1)]
    colnames(py_focal_tokeep_renamed)[colnames(py_focal_tokeep_renamed)=="new_cluster_id"]<-"cluster_id"   ##### py_focal_tokeep_renamed is a pyclone file format ###

    
    ## take the focal variant file and make sudo cluster 0 
    var_focal<-variants_aa[variants_aa$subject_id==subject_id[ss],]   ### take focal variant file 
    
    not_inpy<-var_focal[!(var_focal$pos_id%in%py_focal_tokeep_renamed$pos_id),]  ### find muttaions in the variant file which are not in the pyclone to make sudo cluster
    not_inpy<-not_inpy[!duplicated(not_inpy$pos_id),]

    ### from here
    #not_inpy<-not_inpy[1:10,]

    not_inpy$cluster_id<-0
    not_inpy$cellular_prevalence<-1
    not_inpy$cellular_prevalence_std<-0
    not_inpy$cluster_assignment_prob<-1
    cluster_0<-unique(not_inpy$pos_id)

    out_res_pos<-NULL
    for(rr in 1:length(cluster_0)){   
     #for(rr in 1:5){    ### add fake cluster 0 muttaions and rbind with pyclone  (keep only 5 cluetr 0 muttaions)

        sudo_positions_focal<-not_inpy[not_inpy$pos_id==cluster_0[rr],]
        #sudo_positions_focal<-sudo_positions_focal[1,]
        sudo_positions_focal<-sudo_positions_focal[rep(seq_len(nrow(sudo_positions_focal)), each = length(unique(py_focal_tokeep_renamed$sample_id))), ]   ## repeat each row n= number of samples time 
        sudo_positions_focal<-sudo_positions_focal[,c("Chromosome","Start","Ref","Alt","subject_id","pos_id", "cluster_id", "cellular_prevalence", "cellular_prevalence_std","cluster_assignment_prob","sample_id")]
        
        a<-data.frame("a"=unique(py_focal$sample_rt))
        
        
        sudo_positions_focal<-cbind(sudo_positions_focal,a)
        sudo_positions_focal<-sudo_positions_focal[,c("Chromosome","Start","Ref","Alt","cluster_id","cellular_prevalence","cellular_prevalence_std","cluster_assignment_prob","subject_id","pos_id","a")]

        colnames(sudo_positions_focal)<-c("Chrom","Pos","Ref", "Alt","cluster_id","cellular_prevalence", "cellular_prevalence_std" ,"cluster_assignment_prob" , "subject_id","pos_id","sample_rt")
        out_res_pos<-rbind(sudo_positions_focal,out_res_pos)  
       
      
    }
    py_focal_tokeep_renamed<-py_focal_tokeep_renamed[,c("Chrom" , "Pos" , "Ref", "Alt","cluster_id",  "cellular_prevalence" , "cellular_prevalence_std", "cluster_assignment_prob" ,"subject_id", "pos_id","sample_rt" )]
     
     sudo_and_other<-rbind(out_res_pos,py_focal_tokeep_renamed)
     out_res_var<-rbind(sudo_and_other,out_res_var)
     #out_res_var$cluster_id<-(out_res_var$cluster_id)+1     ## temporary
}


write.table(out_res_var, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/out_res_var_full_pre_post_noRT_pyclone_format_with_cluster0.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


##########################################
##### load phi files from pairtree #######
##########################################
## assign samples in pyclone to their pairtree phi value


out_res_all_phi<-NULL
#for (tt in 1:length(file_trees)){

     for (tt in 1:length(subject_id)){
     
     
     focal_file_tree<-file_trees[tt] 
     json_focal <- fromJSON(paste(readLines(focal_file_tree), collapse=""))  ## parse Json tree file into R

     ## adjust samples name to match with variant file
     subject<-gsub(".tree.json", "",sub(".*/", "", focal_file_tree)) 
     subjecta<-gsub("-.*$","",subject)
     

     #subject<-gsub("-.*$","",subject)

     sample_id<-json_focal$samples
     sample_id<-gsub("R_","",sample_id)


     samples_right<-paste(subjecta,sample_id, sep = "_")
     samples_right

     tree_focal<-data.frame("clone"=0:(length(json_focal$phi)-1))

     out_res_phi<-NULL
     for(jj in 1:nrow(tree_focal)){

          focal_tree_phi<-tree_focal[jj,]
          focal_phi<-data.frame("pairtree_cff"=json_focal$phi[[jj]])

          focal_phi$cluster_id<-focal_tree_phi
          focal_phi$subject_id<-subject
          focal_phi_sample<-cbind(focal_phi,samples_right)
          #focal_phi_sample$pairtree_cff_percent<-(focal_phi_sample$pairtree_cff)*100
        
          out_res_phi<-rbind(focal_phi_sample,out_res_phi)
     }
 out_res_all_phi<-rbind(out_res_phi,out_res_all_phi)     

}    

out_res_all_phi$clust_sample_id<-paste(out_res_all_phi$samples_right , out_res_all_phi$cluster_id , sep = "_")
#    pairtree_cff cluster_id          sample pairtree_cff_percent
#1    6.629058e-01          9 TB9051_1_postRT         6.629058e+01
#2    7.199886e-01          9 TB9051_11_preRT         7.199886e+01
#3    7.412406e-01          9 TB9051_12_preRT         7.412406e+01
out_res_all_phi$subject_id<-gsub("-","_",out_res_all_phi$subject_id)
write.table(out_res_all_phi, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/out_res_full_pre_post_noRT_phi.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##########################################
##########################################
### make clonevol variant file fomat #####
##########################################
##########################################
subject_id<-paired_samples

#out_res_clon<-NULL
for (zz in 1:length(subject_id)){

    out_res_all_phi_focal<-out_res_all_phi[out_res_all_phi$subject_id==subject_id[zz],]  ## pairtree ccf's
    py_focal<-out_res_var[out_res_var$subject_id==subject_id[zz],]    ### out_res_var is the file with sudo-cluster 0 
    py_focal$clust_sample_id<-paste(py_focal$sample_rt , py_focal$cluster_id , sep = "_")
    
    
    #mutaion_RT_focal<-merge(py_focal,meta, by.x = "sample_id", by.y = "sampleid")
    #mutaion_RT_focal$RT_status<-ifelse(mutaion_RT_focal$sequenceofsamplevRT== "beforeRT", "preRT",
    #                    ifelse(mutaion_RT_focal$sequenceofsamplevRT== "nopreopRT", "preRT",
    #                    ifelse(mutaion_RT_focal$sequenceofsamplevRT== "afterRT", "postRT",
    #                    "-")))
    #mutaion_RT_focal$identifier<-paste(mutaion_RT_focal$unique_sample_id.x,mutaion_RT_focal$RT_status, sep = "_")
    #mutaion_RT_focal$code<-paste(mutaion_RT_focal$sample_id,mutaion_RT_focal$RT_status, sep = "_")
    #mutaion_RT_focal$pos_identifier<-gsub("chr","",mutaion_RT_focal$pos_id)

    #pos_identifier<-unique(mutaion_RT_focal$pos_identifier)
    #mutaion_RT_focal_good_col<-mutaion_RT_focal[,c("code","cluster_id","cellular_prevalence","pos_identifier","sample_id","clust_sample_id")]
    
     
     mutaion_RT_focal_good<-merge(py_focal,out_res_all_phi_focal, by.x = "clust_sample_id", by.y = "clust_sample_id")[,c("subject_id.x","cluster_id.x","pairtree_cff","sample_rt","pos_id")]
     pos_identifier<-unique(mutaion_RT_focal_good$pos_id)

     
     out_res_phy<-NULL
     for(ii in 1:length(pos_identifier)) {

          #for(ii in 1:120) {
     

    focal_pos<-mutaion_RT_focal_good[mutaion_RT_focal_good$pos_id== pos_identifier[ii],]
    focal_pos
     
     if(length(unique(focal_pos$out_res_all_phi_focal))< length(unique(out_res_all_phi_focal$sample_rt))){

         find_not_included<-setdiff(unique(out_res_all_phi_focal$sample_rt),focal_pos$sample_rt)  ### find sample doesnot have the fake position
         focal_pos<-rbind(focal_pos, focal_pos[rep(1,length(find_not_included)), ])   ## add rows based on the number of missing samples
         focal_pos$sample_rt[duplicated(focal_pos$sample_rt)] <- find_not_included

    }

    focal_pos$pairtree_cff<-(focal_pos$pairtree_cff)*100
    #focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)/2


    #focal_pos<-out_res_mutaion_RT_focal_good[out_res_mutaion_RT_focal_good$pos_identifier== pos_identifier[ii],]
    #focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)*100

    aa<-data.frame(t(focal_pos))
    colnames(aa)<-aa["sample_rt",]
    bb<-aa["pairtree_cff",]
    rownames(bb)<-NULL
    colnames(bb)<-paste(colnames(bb),"ccf", sep = ".")
    bb[] <- lapply(bb, function(x) as.numeric(as.character(x)))
        
    vaf<-bb/2
    colnames(vaf)<-gsub("ccf","vaf",colnames(vaf))
    
    vaf2<-vaf
    colnames(vaf2)<-gsub(".vaf","",colnames(vaf2))
    
    both<-cbind(vaf,bb,vaf2)
    both$cluster<-unique(focal_pos$cluster_id.x)
    both$position<-unique(focal_pos$pos_identifier)
    variant_focal<-variants_aa[variants_aa$subject_id ==subject_id[zz] & variants_aa$pos_id==pos_identifier[ii],c("gene","is.driver")]
    if(nrow(variant_focal)==0){

     variant_focal<-data.frame("gene"= "-","is.driver" = "FALSE")
    }
    
    both_variant_focal<-cbind(variant_focal,both)
    both_variant_focal$cluster<-(both_variant_focal$cluster)+1


    out_res_phy<-rbind(both_variant_focal,out_res_phy)
    #write.table(file = "out_res_phy", file = paste("~/Desktop/out_res_variants_clonevol_format_",subject_id[zz],".txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
       
       } ## pos_identifier loop
     #out_res_phy$cluster<-(out_res_phy$cluster)+1
     write.table(out_res_phy, file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/variants.",subject_id[zz],".tsv", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

} ##


 

#####################
#####################
### make trees #####
#####################
#####################


file_trees<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/trees", pattern = "*.tsv", full.names = TRUE)
file_vars<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/variants", pattern = "*.tsv", full.names = TRUE)


for (ii in 1:length(file_trees)){

      y = import.tree(file_trees[ii], file_vars[ii])
      y = convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')


      y <- transfer.events.to.consensus.trees(y,
      y$variants[y$variants$is.driver,],
      cluster.col.name = 'cluster',
      event.col.name = 'gene')

      tree_title<-gsub("/Users/shsoudi/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/trees/", "",gsub(".tsv","",file_trees[ii]))


      pdf(file = paste("/Users/shsoudi/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/pairtree_to_clonevol/trees_final_March08/batch2_full/output/clonevol_trees/clonevol_tree_",tree_title,".pdf", sep = ""), width = 7, height = 7)
      plot.all.trees.clone.as.branch(y, branch.width = 0.3,
      node.size = 2.5, node.label.size = 0.5,branch.texts.size = 6)
      dev.off()

}



############
#### END ###
