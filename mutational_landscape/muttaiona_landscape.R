
############################################################
##### process mutations and run dNdScv fpr Sarcoma samples (all samples) #####
#rm(list = ls())
library("tidyverse")
library("dndscv")
library("plyr")
library("VennDiagram")
library("maftools")

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")

## selected samples
#meta<-read.table(file = "selected_samples.txt", header = FALSE, sep= "\t")
#colnames(meta)<-c("sample_id","purity","ploidy","RT_status","RT_status2","RT_staus_code")
#colnames(meta)<-c("sample_id","RT_status","RT_staus_code")

## all samples
meta<-read.table(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)


### load and read all maf files
filenames_CN <- list.files("filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=FALSE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange")]
     })

shared_clonal_mutations <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix<-shared_clonal_mutations%>%drop_na(sample_id)
shared_clonal_mutations_fix$Chromosome<-gsub("chr","",shared_clonal_mutations_fix$Chromosome)

#mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange","RT_status","RT_staus_code")]
mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")
mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]


mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

mutaion_RT$RTtretment<-ifelse(mutaion_RT$RT_status== "preRT", "naive",
                        ifelse(mutaion_RT$RT_status== "noRT", "naive",
                        ifelse(mutaion_RT$RT_status== "postRT", "treatment",
                        "-")))


#cols <- c("Impact","Coding")
#mutaion_RT[cols] <- lapply(mutaion_RT[cols], as.character)                   
#mutaion_RT$Impact[mutaion_RT$Coding== "splicing"] <- "Splice sites"
#colnames(mutaion_RT)<-c("sampleID","chr","pos","ref","mut","RT_code","Gene","Impact","RTtretment","RT_staus_code")

######################################
######## plot base statistics ########
######### all muttaions ##############

mutaion_RT_pre<-mutaion_RT[mutaion_RT$RT_code=="T1",]
mutaion_RT_pre$pos_id<-paste(mutaion_RT_pre$Chromosome,mutaion_RT_pre$Start, sep = "_")
mutaion_RT_pre_mut<-unique(mutaion_RT_pre$pos_id)


mutaion_RT_post<-mutaion_RT[mutaion_RT$RT_code=="T2",]
mutaion_RT_post$pos_id<-paste(mutaion_RT_post$Chromosome,mutaion_RT_post$Start, sep = "_")
mutaion_RT_post_mut<-unique(mutaion_RT_post$pos_id)


mutaion_RT_noRT<-mutaion_RT[mutaion_RT$RT_code=="T3",]
mutaion_RT_noRT$pos_id<-paste(mutaion_RT_noRT$Chromosome,mutaion_RT_noRT$Start, sep = "_")
mutaion_RT_noRT_mut<-unique(mutaion_RT_noRT$pos_id)

### find mutations shared by all three groups
mutaion_RT_pre_post_intersect<-merge(mutaion_RT_pre, mutaion_RT_post, by.x ="pos_id" , by.y = "pos_id")
all_intersect<-merge(mutaion_RT_pre_post_intersect, mutaion_RT_noRT, by.x ="pos_id" , by.y = "pos_id")  ####


## mutations which are unique to preRT ##
RT_pre_mut_unique<-mutaion_RT_pre[!(mutaion_RT_pre$pos_id%in%all_intersect$pos_id),]
RT_pre_mut_unique$plot_code<-rep("T1_preRT")
count(RT_pre_mut_unique$sample_id)
mean_pre<-data.frame("mean_mutations"=mean(count(RT_pre_mut_unique$sample_id)$freq))
#sd(count(RT_pre_mut_unique$sample_id)$freq)
#common_pre<-mutaion_RT_pre[mutaion_RT_pre$pos_id=="1_144816621",]
#common_pre$plot_code<-rep("common_threeway")



RT_post_mut_unique<-mutaion_RT_post[!(mutaion_RT_post$pos_id%in%all_intersect$pos_id),]
RT_post_mut_unique$plot_code<-rep("T2_preRT")
mean_post<-data.frame("mean_mutations"=mean(count(RT_post_mut_unique$sample_id)$freq))
#sd(count(RT_post_mut_unique$sample_id)$freq)
#common_post<-mutaion_RT_post[mutaion_RT_post$pos_id=="1_144816621",]
#common_post$plot_code<-rep("common_threeway")


RT_nort_mut_unique<-mutaion_RT_noRT[!(mutaion_RT_noRT$pos_id%in%all_intersect$pos_id),]
RT_nort_mut_unique$plot_code<-rep("T3_noRT")
mean_noRT<-data.frame("mean_mutations"=mean(count(RT_nort_mut_unique$sample_id)$freq))
#sd(count(RT_nort_mut_unique$sample_id)$freq)
#common_noRT<-mutaion_RT_noRT[mutaion_RT_noRT$pos_id=="1_144816621",]
#common_noRT$plot_code<-rep("common_threeway")


qq<-rbind(mean_pre,mean_post,mean_noRT,1)
qq$rt_code<-c("preRT","postRT","noRT","common")

ggplot(qq) +
  geom_bar( aes(x=rt_code, y=mean_mutations,fill=rt_code), stat="identity",  alpha=0.9)
  #geom_crossbar( aes(x=rt_code, y=mean_mutations, ymin=mean_mutations-sd, ymax=mean_mutations+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)

#all_with_plot_id<-rbind(RT_pre_mut_unique,common_pre,RT_post_mut_unique,common_post,RT_nort_mut_unique,common_noRT)


#mutaion_RT$pos_id<-paste(mutaion_RT$Chromosome,mutaion_RT$Start, sep = "_")
#mutaion_RT_selected<-mutaion_RT[,c("sample_id", "Chromosome","Start","RT_staus_code","pos_id")]
#mut_per_sample<-count(mutaion_RT_selected$sample_id)
#mut_per_sample_withID<-merge(mut_per_sample,meta,by.x = "x", by.y="sample_id")
#summary(mut_per_sample_withID)
#one.way <- aov(freq ~ RT_status, data = mut_per_sample_withID) 

###########################################
####### preRT and postRT comparison #######

#mutaion_RT_pre_post_intersect<-merge(mutaion_RT_pre, mutaion_RT_post, by.x ="pos_id" , by.y = "pos_id")
mutaion_RT_pre_post_intersect<-intersect(mutaion_RT_pre$pos_id, mutaion_RT_post$pos_id)

## mutations which are unique to preRT ##
RT_pre_two_way_uniq<-mutaion_RT_pre[!(mutaion_RT_pre$pos_id%in%mutaion_RT_pre_post_intersect),]
RT_pre_two_way_uniq$plot_code<-rep("T1_preRT")
#count(RT_pre_two_way_uniq$sample_id)
mean_pre_uniq<-data.frame("mean_mutations"=mean(count(RT_pre_two_way_uniq$sample_id)$freq))
count_uniq_pre<-count(RT_pre_two_way_uniq$sample_id)
count_uniq_pre$status<-"preRT"

## mutations which are unique to postRT ##
RT_post_two_way_uniq<-mutaion_RT_post[!(mutaion_RT_post$pos_id%in%mutaion_RT_pre_post_intersect),]
RT_post_two_way_uniq$plot_code<-rep("T2_preRT")
#count(RT_post_two_way_uniq$sample_id)
mean_post_uniq<-data.frame("mean_mutations"=mean(count(RT_post_two_way_uniq$sample_id)$freq))
count_uniq_post<-count(RT_post_two_way_uniq$sample_id)
count_uniq_post$status<-"postRT"


## shared
RT_pre_two_way_shared<-mutaion_RT_pre[(mutaion_RT_pre$pos_id%in%mutaion_RT_pre_post_intersect),]
RT_pre_two_way_shared$plot_code<-rep("shared")
count(RT_pre_two_way_shared$sample_id)
mean_pre_shared<-data.frame("mean_mutations"=mean(count(RT_pre_two_way_shared$sample_id)$freq))


RT_post_two_way_shared<-mutaion_RT_post[(mutaion_RT_post$pos_id%in%mutaion_RT_pre_post_intersect),]
RT_post_two_way_shared$plot_code<-rep("shared")
count(RT_post_two_way_shared$sample_id)
mean_post_shared<-data.frame("mean_mutations"=mean(count(RT_post_two_way_shared$sample_id)$freq))

shared_both<-rbind(RT_pre_two_way_shared,RT_post_two_way_shared)
count(shared_both$sample_id)
mean_both<-data.frame("mean_mutations"=mean(count(shared_both$sample_id)$freq))
count_common<-count(shared_both$sample_id)
count_common$status<-"common"

qq<-rbind(mean_pre_uniq,mean_post_uniq,mean_both)
qq$rt_code<-c("preRT","postRT","common")

ggplot(qq) +
  geom_bar( aes(x=rt_code, y=mean_mutations,fill=rt_code), stat="identity",  alpha=0.9)
  #geom_crossbar( aes(x=rt_code, y=mean_mutations, ymin=mean_mutations-sd, ymax=mean_mutations+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)

all_counts<-rbind(count_uniq_pre,count_uniq_post,count_common)
bp <- ggplot(all_counts, aes(x=status, y=freq, fill=status)) + 
  geom_boxplot()+
  labs(title="Plot of length  per dose",x="Dose (mg)", y = "Length")
bp + theme_classic()



p <- ggplot(all_counts, aes(x=status, y=freq, fill = status)) + 
  geom_boxplot()

p + geom_jitter(shape=16, position=position_jitter(0.2))
##########################################################
##### shared and unique muttaions in pair samples #######
##### samples that have both preRT and postRT samples ####

samples_both<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167")
mutaion_paired<-mutaion_RT[mutaion_RT$unique_sample_id%in%samples_both,]

mutaion_paired$pos_id<-paste(mutaion_paired$Chromosome,mutaion_paired$Start, sep = "_")

out_res<-NULL
for (i in 1:length(samples_both)){

    mutaion_paired_focal<- mutaion_paired[mutaion_paired$unique_sample_id==samples_both[i],]
    mutaion_paired_focal_pre<-mutaion_paired_focal[mutaion_paired_focal$RT_code=="T1",]
    mutaion_paired_focal_post<-mutaion_paired_focal[mutaion_paired_focal$RT_code=="T2",]
#muttaions_shared_pre_post<-merge(mutaion_paired_focal_pre,mutaion_paired_focal_post, by.x = "pos_id",by.y ="pos_id")

    common <- intersect(mutaion_paired_focal_pre$pos_id, mutaion_paired_focal_post$pos_id)  
    pre_uniq<-mutaion_paired_focal_pre[!(mutaion_paired_focal_pre$pos_id%in%common),]
    post_uniq<-mutaion_paired_focal_post[!(mutaion_paired_focal_post$pos_id%in%common),]

    final<-data.frame("common"=length(common))
    final$pre<-length(unique(pre_uniq$pos_id))
    final$post<-length(unique(post_uniq$pos_id))
    final_t<-data.frame(t(final))
    final_t$sharing_status<-rownames(final_t)
    rownames(final_t)<-NULL
    final_t$sample_id<-rep(samples_both[i])
    out_res<-rbind(final_t,out_res)
}

colnames(out_res)<-c("count","sharing_status","sample_id")


p <- ggplot(data = out_res, aes(x = sample_id, y = count)) +
  geom_col(aes(fill = sharing_status), width = 0.7)+
  geom_text(aes( label = count, group =sharing_status), color = "white")
p


ggplot(out_res, aes(x = sample_id, y = count, fill = sharing_status)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#DADAEB", "#9E9AC8", "#6A51A3"))



###########################################
######## only nonsynonymous mutations #####
#### subset ONLy nonsyn mutations #########

mut_nonsyn<-mutaion_RT[(mutaion_RT$Impact=="nonsynonymous SNV" | mutaion_RT$Impact=="stopgain" |mutaion_RT$Impact=="stoploss"),]
mut_nonsyn$pos_id<-paste(mut_nonsyn$Chromosome,mut_nonsyn$Start, sep = "_")
rt_status_code<-c("T1","T2","T3")

pre<-mut_nonsyn[mut_nonsyn$RT_code=="T1",]   ### preRT_nonsyn
pre$pos_id<-paste(pre$Chromosome,pre$Start, sep = "_")
pre_mut<-unique(pre$pos_id)


post<-mut_nonsyn[mut_nonsyn$RT_code=="T2",]   ### postRT_nonsyn
post$pos_id<-paste(post$Chromosome,post$Start, sep = "_")
post_mut<-unique(post$pos_id)


nort<-mut_nonsyn[(mut_nonsyn$RT_code=="T3" | mut_nonsyn$RT_code=="T4"),]    ### noRT_nonsyn
nort$pos_id<-paste(nort$Chromosome,nort$Start, sep = "_")
nort_mut<-unique(nort$pos_id)


pre_post_intersect<-merge(pre, post, by.x ="pos_id" , by.y = "pos_id")
all_intersect_nonsyn<-merge(pre_post_intersect, nort, by.x ="pos_id" , by.y = "pos_id")  ####


pre_nonsyn_unique<-pre[!(pre$pos_id%in%all_intersect_nonsyn$pos_id),]
#RT_pre_mut_unique$plot_code<-rep("T1_preRT")
mean_pre<-data.frame("mean_mutations"=mean(count(pre_nonsyn_unique$sample_id)$freq))


post_nonsyn_unique<-post[!(post$pos_id%in%all_intersect_nonsyn$pos_id),]
#RT_pre_mut_unique$plot_code<-rep("T1_preRT")
mean_post<-data.frame("mean_mutations"=mean(count(post_nonsyn_unique$sample_id)$freq))


noRT_nonsyn_unique<-nort[!(nort$pos_id%in%all_intersect_nonsyn$pos_id),]
#RT_pre_mut_unique$plot_code<-rep("T1_preRT")
mean_nort<-data.frame("mean_mutations"=mean(count(noRT_nonsyn_unique$sample_id)$freq))

ww<-rbind(mean_pre,mean_post,mean_nort,0)
ww$rt_code<-c("preRT","postRT","noRT","common")

ggplot(ww) +
  geom_bar( aes(x=rt_code, y=mean_mutations,fill=rt_code), stat="identity",  alpha=0.9)


##########################################################
#### Shared and unique beween pre and post treatment #####

pre_post_intersect<-intersect(pre, post)
## mutations which are unique to preRT ##
RT_pre_two_way_uniq<-pre[!(pre$pos_id%in%pre_post_intersect$pos_id),]
RT_pre_two_way_uniq$plot_code<-rep("T1_preRT")
count(RT_pre_two_way_uniq$sample_id)
mean_pre_uniq<-data.frame("mean_mutations"=mean(count(RT_pre_two_way_uniq$sample_id)$freq))


RT_post_two_way_uniq<-post[!(post$pos_id%in%pre_post_intersect$pos_id),]
RT_post_two_way_uniq$plot_code<-rep("T2_preRT")
count(RT_post_two_way_uniq$sample_id)
mean_post_uniq<-data.frame("mean_mutations"=mean(count(RT_post_two_way_uniq$sample_id)$freq))

## shared
RT_pre_two_way_shared<-pre[(pre$pos_id%in%pre_post_intersect$pos_id),]
RT_pre_two_way_shared$plot_code<-rep("shared")
count(RT_pre_two_way_shared$sample_id)
mean_pre_shared<-data.frame("mean_mutations"=mean(count(RT_pre_two_way_shared$sample_id)$freq))


RT_post_two_way_shared<-post[(post$pos_id%in%pre_post_intersect$pos_id),]
RT_post_two_way_shared$plot_code<-rep("shared")
count(RT_post_two_way_shared$sample_id)
mean_post_shared<-data.frame("mean_mutations"=mean(count(RT_post_two_way_shared$sample_id)$freq))

shared_both<-rbind(RT_pre_two_way_shared,RT_post_two_way_shared)
count(shared_both$sample_id)
mean_both<-data.frame("mean_mutations"=mean(count(shared_both$sample_id)$freq))


qq<-rbind(mean_pre_shared,mean_post_shared,mean_both)
qq$rt_code<-c("preRT","postRT","common")

ggplot(qq) +
  geom_bar( aes(x=rt_code, y=mean_mutations,fill=rt_code), stat="identity",  alpha=0.9)
  #geom_crossbar( aes(x=rt_code, y=mean_mutations, ymin=mean_mutations-sd, ymax=mean_mutations+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)

#################################################
#################################################
######## different baepair type #################
## find basechanges "nonsynonymous only"
unique(mut_nonsyn$basechange)
#G|A T|A C|T G|T C|A G|C A|C A|T T|C C|G T|G A|G

u_basechange<-c("C|A","C|G","C|T","T|A","T|C","T|G")
rt_status<-unique(mut_nonsyn$RT_status)


mut_nonsyn$base[mut_nonsyn$basechange== "G|A"] <- "C|T"
mut_nonsyn$base[mut_nonsyn$basechange== "G|T"] <- "C|A"
mut_nonsyn$base[mut_nonsyn$basechange== "G|C"] <- "C|G"
mut_nonsyn$base[mut_nonsyn$basechange== "A|C"] <- "T|G"
mut_nonsyn$base[mut_nonsyn$basechange== "A|T"] <- "T|A"
mut_nonsyn$base[mut_nonsyn$basechange== "A|G"] <- "T|C"



mutaion_RT$Variant_Classification[mutaion_RT$Impact== "stopgain"] <- "Nonsense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "nonsynonymous SNV"] <- "Missense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "frameshift insertion"] <- "Frame_Shift_Ins"



our_res<-NULL
for (i in 1:length(u_basechange)){

  focal<-mut_nonsyn[mut_nonsyn$basechange==u_basechange[i],]

  for (j in 1:length(rt_status)){

    focal_rt<-focal[focal$RT_status==rt_status[j],]
    mut_counts<-count(focal_rt$sample_id)
    mut_counts$mu_type<-u_basechange[i]
    mut_counts$rt_status<-rt_status[j]
    our_res<-rbind(mut_counts,our_res)
  }
  
}


ggplot(data = our_res, aes(x = mu_type, y = freq, fill = rt_status, color = rt_status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))


### Transversion versus transition muttaions 
#https://www.differencebetween.com/difference-between-transition-and-vs-transversion/


#########################################################
## limit to pre and postRT samples (filtered out noRT) ##
mut_nonsyn_both_pre_post<-mut_nonsyn[(mut_nonsyn$RT_status=="preRT" | mut_nonsyn$RT_status=="postRT"),]
rt_status<-unique(mut_nonsyn_both_pre_post$RT_status)

our_res_pre_post<-NULL
for (i in 1:length(u_basechange)){

  focal<-mut_nonsyn_both_pre_post[mut_nonsyn_both_pre_post$basechange==u_basechange[i],]
  
  for (j in 1:length(rt_status)){

    focal_rt<-focal[focal$RT_status==rt_status[j],]
    mut_counts<-count(focal_rt$sample_id)
    mut_counts$mu_type<-u_basechange[i]
    mut_counts$rt_status<-rt_status[j]
    our_res_pre_post<-rbind(mut_counts,our_res_pre_post)
  }
  
}

res_TC<- our_res_pre_post[our_res_pre_post$mu_type=="T|G",] 
res_TC_test <- wilcox.test(freq ~ rt_status, data = res_TC)
res_TC_test
boxplot(freq~rt_status, data = res_TC)

ggplot(data = our_res_pre_post, aes(x = mu_type, y = freq, fill = rt_status, color = rt_status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))


#################################################
######## different baepair type #################
## find basechanges "nonsynonymous only"
unique(mutaion_RT$basechange)
#G|A T|A C|T G|T C|A G|C A|C A|T T|C C|G T|G A|G

u_basechange<-c("C|A","C|G","C|T","T|A","T|C","T|G")
rt_status<-unique(mut_nonsyn$RT_status)


mut_nonsyn$base[mut_nonsyn$basechange== "G|A"] <- "C|T"
mut_nonsyn$base[mut_nonsyn$basechange== "G|T"] <- "C|A"
mut_nonsyn$base[mut_nonsyn$basechange== "G|C"] <- "C|G"
mut_nonsyn$base[mut_nonsyn$basechange== "A|C"] <- "T|G"
mut_nonsyn$base[mut_nonsyn$basechange== "A|T"] <- "T|A"
mut_nonsyn$base[mut_nonsyn$basechange== "A|G"] <- "T|C"



mutaion_RT$Variant_Classification[mutaion_RT$Impact== "stopgain"] <- "Nonsense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "nonsynonymous SNV"] <- "Missense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "frameshift insertion"] <- "Frame_Shift_Ins"



our_res<-NULL
for (i in 1:length(u_basechange)){

  focal<-mut_nonsyn[mut_nonsyn$basechange==u_basechange[i],]

  for (j in 1:length(rt_status)){

    focal_rt<-focal[focal$RT_status==rt_status[j],]
    mut_counts<-count(focal_rt$sample_id)
    mut_counts$mu_type<-u_basechange[i]
    mut_counts$rt_status<-rt_status[j]
    our_res<-rbind(mut_counts,our_res)
  }
  
}


ggplot(data = our_res, aes(x = mu_type, y = freq, fill = rt_status, color = rt_status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))


## limit to pre and postRT samples (filtered out noRT)
mut_nonsyn_both_pre_post<-mut_nonsyn[(mut_nonsyn$RT_status=="preRT" | mut_nonsyn$RT_status=="postRT"),]
rt_status<-unique(mut_nonsyn_both_pre_post$RT_status)

our_res_pre_post<-NULL
for (i in 1:length(u_basechange)){

  focal<-mut_nonsyn_both_pre_post[mut_nonsyn_both_pre_post$basechange==u_basechange[i],]
  
  for (j in 1:length(rt_status)){

    focal_rt<-focal[focal$RT_status==rt_status[j],]
    mut_counts<-count(focal_rt$sample_id)
    mut_counts$mu_type<-u_basechange[i]
    mut_counts$rt_status<-rt_status[j]
    our_res_pre_post<-rbind(mut_counts,our_res_pre_post)
  }
  
}

res_TC<- our_res_pre_post[our_res_pre_post$mu_type=="T|G",] 
res_TC_test <- wilcox.test(freq ~ rt_status, data = res_TC)
res_TC_test
boxplot(freq~rt_status, data = res_TC)

ggplot(data = our_res_pre_post, aes(x = mu_type, y = freq, fill = rt_status, color = rt_status)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))















############################################################
#######Find recurrent mutations by running dNdSCV ##########
############################################################
#mmm<-mm[,c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2" ,"Tumor_Sample_Barcode")]       
#unique(mmm$Variant_Classification)
# [1] Intron                 IGR                    3'UTR                  5'Flank               
# [5] 3'Flank                Missense_Mutation      Nonsense_Mutation      RNA                   
# [9] Silent                 5'UTR                  Nonstop_Mutation       Splice_Region         
#[13] Splice_Site            Frame_Shift_Del        Frame_Shift_Ins        In_Frame_Del          
#[17] In_Frame_Ins           Translation_Start_Site

meta<-read.delim(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1,5,6)]
colnames(meta)<-c("sample_id","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)


### load and read all maf files
filenames_CN <- list.files("filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.csv(x, header=TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange")]
     })

shared_clonal_mutations <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix<-shared_clonal_mutations%>%drop_na(sample_id)
shared_clonal_mutations_fix$Chromosome<-gsub("chr","",shared_clonal_mutations_fix$Chromosome)

#mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange","RT_status","RT_staus_code")]
mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")
mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]


mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

#mutaion_RT$RTtretment<-ifelse(mutaion_RT$RT_status== "preRT", "naive",
#                        ifelse(mutaion_RT$RT_status== "noRT", "naive",
#                        ifelse(mutaion_RT$RT_status== "postRT", "treatment",
#                       "-")))


### adjust maf file format               
mutaion_RT$Variant_Classification[mutaion_RT$Coding== "splicing"] <- "Splice_Site"
mutaion_RT$Variant_Classification[mutaion_RT$Coding== "UTR5"] <- "5'UTR"
mutaion_RT$Variant_Classification[mutaion_RT$Coding== "UTR3"] <- "3'UTR"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "synonymous SNV"] <- "Silent"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "stopgain"] <- "Nonsense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "nonsynonymous SNV"] <- "Missense_Mutation"
mutaion_RT$Variant_Classification[mutaion_RT$Impact== "frameshift insertion"] <- "Frame_Shift_Ins"

mutaion_RT<-mutaion_RT[!(mutaion_RT$Impact=="unknown"),]
mutaion_RT_good<-mutaion_RT%>%drop_na(Variant_Classification)


#mutaion_RT[cols] <- lapply(mutaion_RT[cols], as.character)

mutaion_RT_good<-mutaion_RT_good[,c("sample_id","Chromosome", "Start","End","Ref","Alt", "Gene","RT_status", "RT_code","Variant_Classification")]
colnames(mutaion_RT_good)<-c("Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Hugo_Symbol","RT_status", "RT_code","Variant_Classification")
mutaion_RT_good<-mutaion_RT_good[,c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode","RT_status", "RT_code")]

mutaion_RT_good$Variant_Type<-rep("SNP") #### canbe used as maf file 
mutaion_RT_good$Variant_Type[mutaion_RT_good$Variant_Classification=="Frame_Shift_Ins"]<-"INS"

#### save prepared maf like file ####
write.table(mutaion_RT_good,file = "all_samples_updated_converted_like.maf", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
### Run dnds

mutaion_for_dnds<-mutaion_RT_good[,c("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2")]
colnames(mutaion_for_dnds)<-c("sampleID","chr","pos","ref","mut")

  dndsout = dndscv(mutaion_for_dnds)
  sel_cv = dndsout$sel_cv
  signif_genes = sel_cv[sel_cv$qallsubs_cv<0.01, c("gene_name","qallsubs_cv")]
  write.table(signif_genes ,file ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_sel_cv_fdr0.01_all_samples.table", ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
   

  globaldnds_con<-(dndsout$globaldnds)
  write.table(globaldnds_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_globaldnds_fdr0.01_all_samples.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
  
  annotmuts_con<-(dndsout$annotmuts)
  write.table(annotmuts_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_annotmuts_fdr0.01_all_samples.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE) 
  
  genemuts_con<-(dndsout$genemuts)
  write.table(genemuts_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_genemuts_fdr0.01_all_sample.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
  
  mle_submodel_con<-(dndsout$mle_submodel)
  write.table(mle_submodel_con,file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_mle_submodel_fdr0.01_all_sample.table" ,col.names = TRUE, row.names = FALSE, sep = "\t" , quote  = FALSE)
  


######################## MafTools ##############
#################################################
#### analyze the data with maftools #############
rm(list = ls())
library("maftools")
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")
mutaion_RT_good<-read.delim(file = "all_samples_updated_converted_like.maf", header = TRUE)
#mutaion_RT_good<-mutaion_RT_good[!(mutaion_RT_good$Reference_Allele=="-"),]

sar_laml = read.maf(maf = mutaion_RT_good)

#ellis_BD <- shared_clonal_mutations_fix[order(shared_clonal_mutations_fix$sampleID,shared_clonal_mutations_fix$chr,shared_clonal_mutations_fix$pos),]
#ind = which(diff(ellis_BD $Start_Position)==1)
#ellis_consecutive<-ellis_BD [unique(sort(c(ind,ind+1))),]
#toDelete <- seq(1, nrow(ellis_consecutive), 2)
#toDelete_df<-ellis_consecutive[toDelete ,]

#Shows sample summry.
getSampleSummary(sar_laml)
#Shows gene summary.
getGeneSummary(sar_laml)
#shows clinical data associated with samples
getClinicalData(sar_laml)
#Shows all fields in MAF
getFields(sar_laml)#
#Writes maf summary to an output file with basename laml.
#write.mafSummary(maf = sar_laml, basename = 'laml')

###Plotting MAF summary
plotmafSummary(maf = sar_laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

### draw oncolplot
sig_genes<-read.table(file ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/dnds/dnds_sel_cv_fdr0.01_all_samples.table", header = TRUE)
sar_aml_genes<-sig_genes[!grepl("RBM25",sig_genes$gene_name),]
sar_aml_genes<-sar_aml_genes[!grepl("MYO3A",sar_aml_genes$gene_name),]
sar_aml_genes<-sar_aml_genes[!grepl("OR10G8",sar_aml_genes$gene_name),]
sar_aml_genes<-sar_aml_genes[!grepl("HERC1",sar_aml_genes$gene_name),]



sar_aml_genes<-sar_aml_genes[1:50,1]
#sar_aml_genes<-sar_aml_genes[-c("RBM25","MYO3A","OR10G8","HERC1")])




oncoplot(maf = sar_laml, genes = sar_aml_genes)

#laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
#plotTiTv(res = laml.titv)


maf_meta<-mutaion_RT_good[,c("Tumor_Sample_Barcode","RT_status","RT_code")]
colnames(maf_meta)<-c("Tumor_Sample_Barcode","RT_status","RT_code")

sar_laml = read.maf(maf = mutaion_RT_good, clinicalData = maf_meta)

oncoplot(maf = sar_laml, genes = sar_aml_genes, clinicalFeatures = 'RT_status',sortByAnnotation = TRUE,fontSize = 0.4, top = 2)

rainfallPlot(maf = sar_laml, detectChangePoints = TRUE, pointSize = 0.4)



#gisticOncoPlot(gistic = laml.gistic, clinicalData = getClinicalData(x = laml), clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10)

#Selected AML driver genes

#Variant allele frequcnies (Right bar plot)
aml_genes_vaf = subsetMaf(maf = laml, genes = aml_genes, fields = "i_TumorVAF_WU", mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
colnames(aml_genes_vaf)[2] = "VAF"
head(aml_genes_vaf)
#>    Hugo_Symbol      VAF
#> 1:       ASXL1 37.11250
#> 2:       CEBPA 22.00235
#> 3:      DNMT3A 43.51556
#> 4:      DNMT3B 37.14000
#> 5:        EZH2 68.88500
#> 6:        FLT3 34.60294

#MutSig results (Right bar plot)
laml.mutsig = system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
laml.mutsig = data.table::fread(input = laml.mutsig)[,.(gene, q)]
laml.mutsig[,q := -log10(q)] #transoform to log10
head(laml.mutsig)
#>      gene        q
#> 1:   FLT3 12.64176
#> 2: DNMT3A 12.64176
#> 3:   NPM1 12.64176
#> 4:   IDH2 12.64176
#> 5:   IDH1 12.64176
#> 6:   TET2 12.64176


oncoplot(
  maf = laml,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig,
  rightBarLims = c(0, 20)
)

