## muttaional signature analysis for merged regions

rm(list = ls())


library(dplyr)
library(reshape2)
library(kableExtra)
library(ggplot2)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)

# Load mutSignatures
library(mutSignatures)
# prep hg19
hg19 <- BSgenome.Hsapiens.UCSC.hg19
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")

### make input file for mutsignatures ###
### load and read all maf files
filenames_CN <- list.files("filtered_variants/selected",pattern="*.filtered.variants.oxomerge.final.txt_selected", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.delim(x, header=FALSE, sep = "\t")
     })

### add sample name as a column
for (i in 1:length(attackStats_CN)){
    attackStats_CN[[i]]<-cbind(attackStats_CN[[i]],filenames_CN[i])
    }
aa <- do.call("rbind", attackStats_CN) 
aa$sample_id<-gsub("filtered_variants/selected/", "",gsub(".filtered.variants.oxomerge.final.txt_selected","",aa[,10]))
shared_clonal_mutations_fix<-aa[-1,c(-10,-8)]
colnames(shared_clonal_mutations_fix)<-c("Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","sample_id")

#shared_clonal_mutations_fix<-shared_clonal_mutations%>%drop_na(sample_id)
shared_clonal_mutations_fix$Chromosome<-gsub("chr","",shared_clonal_mutations_fix$Chromosome)
shared_clonal_mutations_fix$pos_id<-paste(shared_clonal_mutations_fix$Chromosome,shared_clonal_mutations_fix$Start, sep = "_")

## all samples
meta<-read.table(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)



mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")
mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]


mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

mutaion_RT$RTtretment<-ifelse(mutaion_RT$RT_status== "preRT", "naive",
                        ifelse(mutaion_RT$RT_status== "noRT", "naive",
                        ifelse(mutaion_RT$RT_status== "postRT", "treatment",
                        "-")))

mutaion_RT<-mutaion_RT[mutaion_RT$Chromosome!="Start",]
#cols <- c("Impact","Coding")
#mutaion_RT[cols] <- lapply(mutaion_RT[cols], as.character)                   
#mutaion_RT$Impact[mutaion_RT$Coding== "splicing"] <- "Splice sites"
#colnames(mutaion_RT)<-c("sampleID","chr","pos","ref","mut","RT_code","Gene","Impact","RTtretment","RT_staus_code")

mutaion_RT$pos_id<-paste(mutaion_RT$Chromosome,mutaion_RT$Start, sep = "_")
mutaion_RT$identifier<-paste(mutaion_RT$unique_sample_id,mutaion_RT$RT_status, sep = "_")


samples<-unique(mutaion_RT$identifier)

out_res<-NULL
for (i in 1:length(samples)){

    mutaion_RT_focal<-mutaion_RT[mutaion_RT$identifier==samples[i],]
    mutaion_RT_focal_nodup<-mutaion_RT_focal[!(duplicated(mutaion_RT_focal$pos_id)),]
    out_res<-rbind(mutaion_RT_focal_nodup,out_res)  
}

out_res_good<-out_res[,c("Chromosome", "Start","Ref","Alt","identifier")]
colnames(out_res_good)<-c("CHROM","POS","REF","ALT","SAMPLEID")
out_res_good$CHROM<-paste("chr",out_res_good$CHROM,sep = "")

## write input file for mutsignatures ###
write.table(out_res_good, file = "outres_input_mutsignatures_merged_regions_all_samples_tab.table", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



input_sigmutation<-read.table(file = "MutSignatures/outres_input_mutsignatures_merged_regions_all_samples_tab.table", header = TRUE)
### filter out any muttaion that is not a SNV

muts<-c("A","T","C","G")
input_sigmutation_snv<-input_sigmutation[input_sigmutation$REF%in%muts,]
input_sigmutation_snv<-input_sigmutation_snv[input_sigmutation_snv$ALT%in%muts,]
#write.table(input_sigmutation_snv, file = "outres_input_mutsignatures_merged_regions_all_samples_tab.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

########################################################################
########################################################################
########################################################################
#### start with Mutsignatures ####
#### extrcat preRT samples
#y_pre<-input_sigmutation_snv[grepl("noRT",input_sigmutation_snv$SAMPLEID),]
y_pre<-input_sigmutation_snv
#De novo extraction of Mutational Signatures from BLCA samples
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



write.table(y_pre, file = "MutSignatures/outres_compute_mutType_MutSignatures_all_samples_merged_regions.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

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

num.sign <- 4

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
pre.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 500,               # bootstrapping: usually 500-1000
    num_parallelCores = 6)                # total num of cores to use (parallelization)


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
msigPlot(pre.sig, signature = 1, ylim = c(0, 0.10))
msigPlot(pre.sig, signature = 2, ylim = c(0, 0.10))
msigPlot(pre.sig, signature = 3, ylim = c(0, 0.10))
msigPlot(pre.sig, signature = 4, ylim = c(0, 0.10))



msigPlot(pre.exp) + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))
xprt <- coerceObj(x = pre.exp, to = "data.frame") 
#xprt <- coerceObj(x = pre.exp, to = "data.frame") 

write.table(xprt, file = "MutSignatures/out_put_mutsignatures_xprt_denove_signature_count_All_RT.txt")


# Retrieve COSMIC signatures from online repo, and then subset
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
cosmx<-cosmix[c(1, 2, 5,13)]

# match OVcar and COSMIC signatures
mSign.sar <- matchSignatures(mutSign = pre.sig, reference = cosmix)
print(mSign.sar$plot)


blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = pre.sig)

blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.1x_df <- coerceObj(x = blca.exp.1x, to = "data.frame") 
write.table(blca.exp.1x_df, file = "MutSignatures/out_put_mutsignatures_blca.exp.1x_cosmic_signature_count_All_RT.txt")


blca.exp.2x <- blca.expo2$Results$count.result
blca.exp.2x_df <- coerceObj(x = blca.exp.2x, to = "data.frame") 
write.table(blca.exp.2x_df, file = "MutSignatures/out_put_mutsignatures_blca.exp.2x_denove_signature_count_All_RT.txt")


# Plot exposures
bp1 <- msigPlot(blca.exp.1x, main = "BLCA | COSMIC sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

bp2 <- msigPlot(blca.exp.2x, main = "BLCA | pre-blca sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

# Visualize
grid.arrange(bp1, bp2, ncol = 2)


#################################
############# postRT #############


post<-mutaion_RT[mutaion_RT$RT_status=="postRT",]
y_post<-post[,c("Chromosome","Start","Ref","Alt","sample_id")]
colnames(y_post)<-c("CHROM","POS","REF","ALT","SAMPLEID")


#De novo extraction of Mutational Signatures from BLCA samples
# Attach context
y_post <- attachContext(mutData = y_post,
                   chr_colName = "CHROM",
                   start_colName = "POS",
                   end_colName = "POS",
                   nucl_contextN = 3,
                   BSGenomeDb = hg19)



# Remove mismatches
y_post <- removeMismatchMut(mutData = y_post,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")    
# Compute mutType
y_post <- attachMutType(mutData = y_post,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context") 


#
post.counts <- countMutTypes(mutTable = y_post,
                             mutType_colName = "mutType",
                             sample_colName = "SAMPLEID")
# Mutation Counts
print(post.counts)                  

 mouCancer.assess <- prelimProcessAssess(input = post.counts, approach = "counts")


#sarcoma.counts.data <- getCounts(pre.counts)
x <- getCounts(post.counts)
#ocd <- as.mutation.counts(sarcoma.counts.data)
xx <- as.mutation.counts(x)
# how many signatures should we extract? 

mouCancer.assess <- prelimProcessAssess(input = xx, approach = "counts")

num.sign <- 4

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
post.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 100,               # bootstrapping: usually 500-1000
    num_parallelCores = 6)                # total num of cores to use (parallelization)


# Extract new signatures (de-novo mutations)- may take a while
post.analysis <- 
  decipherMutationalProcesses(input = post.counts,
                              params = post.params)



#Downstream analyses and visualization
## examine the results
# Retrieve signatures (results)
post.sig <- post.analysis$Results$signatures

# Retrieve exposures (results)
post.exp <- post.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(post.sig, signature = 1, ylim = c(0, 0.10))
msigPlot(post.sig, signature = 2, ylim = c(0, 0.10))
msigPlot(post.sig, signature = 3, ylim = c(0, 0.10))
msigPlot(post.sig, signature = 4, ylim = c(0, 0.10))
#msigPlot(pre.sig, signature = 5, ylim = c(0, 0.10))




msigPlot(post.exp) + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))


# Export Signatures as data.frame
xprt <- coerceObj(x = post.sig, to = "data.frame") 
write.table(xprt, file = "")

# Get signatures from data (imported as data.frame) 
# and then convert it to mutSignatures object


# Retrieve COSMIC signatures from online repo, and then subset
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
cosmx<-cosmix[c(1, 2, 5,13)]


# match OVcar and COSMIC signatures
mSign.sar <- matchSignatures(mutSign = post.sig, reference = cosmix)
print(mSign.sar$plot)


#########
#pre.sig -> blcmx
#cosmix

blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = post.sig)




blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.2x <- blca.expo2$Results$count.result

# Plot exposures
bp1 <- msigPlot(blca.exp.1x, main = "BLCA | COSMIC sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

bp2 <- msigPlot(blca.exp.2x, main = "BLCA | pre-blca sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

# Visualize
grid.arrange(bp1, bp2, ncol = 2)




#####



qq<-data.frame(t(blca.exp.2x_df))

aa<-data.frame(qq[,1])
aa$sample_id<-rownames(qq)
aa$sig<-rep("sig.1")
colnames(aa)<-c("count","sample_id","sig")


bb<-data.frame(qq[,2])
bb$sample_id<-rownames(qq)
bb$sig<-rep("sig.2")
colnames(bb)<-c("count","sample_id","sig")


cc<-data.frame(qq[,3])
cc$sample_id<-rownames(qq)
cc$sig<-rep("sig.3")
colnames(cc)<-c("count","sample_id","sig")



dd<-data.frame(qq[,4])
dd$sample_id<-rownames(qq)
dd$sig<-rep("sig.4")
colnames(dd)<-c("count","sample_id","sig")

all<-rbind(aa,bb,cc,dd)

pdf("~/desktop/cexp_counts_annotated.pdf")
plota<-ggplot(all, aes(x = sample_id, y= count,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plota)
dev.off()   



####




vv<-data.frame(t(xprt))

aa<-data.frame(vv[,1])
aa$sample_id<-rownames(vv)
aa$sig<-rep("sig.1")
colnames(aa)<-c("count","sample_id","sig")


bb<-data.frame(vv[,2])
bb$sample_id<-rownames(vv)
bb$sig<-rep("sig.2")
colnames(bb)<-c("count","sample_id","sig")


cc<-data.frame(vv[,3])
cc$sample_id<-rownames(vv)
cc$sig<-rep("sig.3")
colnames(cc)<-c("count","sample_id","sig")



dd<-data.frame(vv[,4])
dd$sample_id<-rownames(vv)
dd$sig<-rep("sig.4")
colnames(dd)<-c("count","sample_id","sig")

all_me<-rbind(aa,bb,cc,dd)

pdf("~/desktop/cexp_counts_annotated.pdf")
plotv<-ggplot(all_me, aes(x = sample_id, y= count,fill = sig)) + 
   geom_bar(stat = "identity")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plotv)
dev.off()   