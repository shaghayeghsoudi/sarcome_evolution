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

### with all samples in the metadata
meta<-read.table(file = "metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
#colnames(meta)<-c("sample_id","purity","ploidy","RT_status","RT_status2","RT_staus_code")
#colnames(meta)<-c("sample_id","RT_status","RT_staus_code")
#meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)


### load and read all maf files
filenames_CN <- list.files("filtered_variants",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.csv(x, header=TRUE, sep = "\t")[,c("sample_id","Chromosome","Start" ,"End","Ref","Alt","Coding","Gene","Impact","V4","V5","basechange")]
     })

shared_clonal_mutations <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix<-shared_clonal_mutations%>%drop_na(sample_id)
mutaion_RT<-merge(shared_clonal_mutations_fix,meta, by.x = "sample_id", by.y = "sample_id")
mutaion_RT<-mutaion_RT[mutaion_RT$Start!="777428",]

mutaion_RT$RT_status<-gsub("noRT-1","noRT",mutaion_RT$RT_status)
mutaion_RT$RT_status<-gsub("noRT-2","noRT",mutaion_RT$RT_status)

#### prep to run signature analysis
### preRT samples
# Data Import and pre-processing
pre<-mutaion_RT[mutaion_RT$RT_status=="preRT",]
y_pre<-pre[,c("Chromosome","Start","Ref","Alt","sample_id")]
colnames(y_pre)<-c("CHROM","POS","REF","ALT","SAMPLEID")

#CHROM      POS REF ALT SAMPLEID
#1 chr10 50863185   G   A SRC125_1
#2  chr7 87170740   T   A SRC125_1
#3 chr12 10370686   G   A SRC125_1
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
    num_totIterations = 100,               # bootstrapping: usually 500-1000
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
#msigPlot(pre.sig, signature = 5, ylim = c(0, 0.10))




msigPlot(pre.exp) + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))


# Export Signatures as data.frame
xprt <- coerceObj(x = pre.sig, to = "data.frame") 
write.table(xprt, file = "")

# Get signatures from data (imported as data.frame) 
# and then convert it to mutSignatures object


# Retrieve COSMIC signatures from online repo, and then subset
cosmix <- getCosmicSignatures() 
#cosmx <- as.mutation.signatures(cosmx)
cosmx<-cosmix[c(1, 2, 5,13)]


# match OVcar and COSMIC signatures
mSign.sar <- matchSignatures(mutSign = pre.sig, reference = cosmix)
print(mSign.sar$plot)


#########
#pre.sig -> blcmx
#cosmix

blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = pre.sig)




blca.exp.1x <- blca.expo1$Results$count.result
blca.exp.2x <- blca.expo2$Results$count.result

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