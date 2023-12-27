## make input "vcf" file for MRT sotware

## Load required libraries
library(tidyverse)

## Example of vcf input file format

##fileformat=VCFv4.1 
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (number of reads in each observed allele)"> 
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth"> 
##FORMAT=<ID=FT,Number=1,Type=String,Description="Variant filters"> 
##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description="Tumour ref count"> 
##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description="Tumour alt count"> 
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
#5	158636	.	C	T	.	.	t_alt_count=4;t_ref_count=10	AD:DP	.,4:14
#5	1758383	.	A	T	.	.	t_alt_count=5;t_ref_count=32	AD:DP	.,5:37


setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs")

## This piece of code converts MAF files (or variant files) into a vcf file format required by the software
## load variant files
filenames_CN <- list.files("variantfiles_test",pattern="*filtered.variants.oxomerge.final.txt", full.names = TRUE)
attackStats_CN <- lapply(filenames_CN,function(x) {
     read.csv(x, header=TRUE, sep = "\t")[,c("sample_id","Chromosome","Start","Ref", "Alt","Gene","AF","TumorReads","TumorDepth")]
     })


shared_clonal_mutations <- do.call("rbind", attackStats_CN) 
shared_clonal_mutations_fix<-shared_clonal_mutations[!rowSums(nchar(as.matrix(shared_clonal_mutations[4:5]))!=1),]


shared_clonal_mutations_fix$Chromosome<-gsub("chr","",shared_clonal_mutations_fix$Chromosome)
shared_clonal_mutations_fix$ID<-rep(".",nrow(shared_clonal_mutations_fix))
shared_clonal_mutations_fix$QUAL<-rep(".",nrow(shared_clonal_mutations_fix))
shared_clonal_mutations_fix$FILTER<-rep(".",nrow(shared_clonal_mutations_fix))
shared_clonal_mutations_fix$FORMAT<-rep("AD:DP",nrow(shared_clonal_mutations_fix))
shared_clonal_mutations_fix$Ref_count<-(shared_clonal_mutations_fix$TumorDepth-shared_clonal_mutations_fix$TumorReads)
shared_clonal_mutations_fix$allele<-paste(".",(paste(shared_clonal_mutations_fix$TumorReads,shared_clonal_mutations_fix$Ref_count, sep =":")),sep=",")

shared_clonal_mutations_fix$Ref_count<-paste("t_ref_count=",shared_clonal_mutations_fix$Ref_count,sep = "")
shared_clonal_mutations_fix$TumorReads<-paste("t_alt_count=",shared_clonal_mutations_fix$TumorReads,sep = "")
shared_clonal_mutations_fix$INFO<-paste(shared_clonal_mutations_fix$TumorReads,shared_clonal_mutations_fix$Ref_count, sep =";" )

muttaion_selected_cols<-shared_clonal_mutations_fix[,c("Chromosome","Start","ID","Ref","Alt","QUAL","FILTER","INFO","FORMAT","allele","sample_id")]
colnames(muttaion_selected_cols)<-c("#CHROM","POS","ID","REF","ALT"	,"QUAL"	,"FILTER",	"INFO",	"FORMAT","SAMPLE", "sample_id")


prettyprint <- function() {
  cat('##fileformat=VCFv4.1',
  '\n##FORMAT=<ID=AD,Number=2,Type=Integer,Description= "Allelic depths (number of reads in each observed allele)" >',
  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth" >',
  '\n##FORMAT=<ID=FT,Number=1,Type=String,Description="Variant filters" >',
  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Tumour ref count" >',
  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Tumour alt counth" > \n'  
      )
      print(muttaion_selected_cols)
 }
 prettyprint()


#########################################
# worked
#prettyprint <- function() {
#     print(df)
#     cat('##FORMAT=<ID=AD,Number=2,Type=Integer,Description= "Allelic depths (number of reads in each observed allele)" >',
#         '\nFORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (number of reads in each observed allele)" >',
#         "\n##%GC=30")
# }
# prettyprint()

#worked

#prettyprint <- function() {
#  cat('##fileformat=VCFv4.1',
#  '\n##FORMAT=<ID=AD,Number=2,Type=Integer,Description= "Allelic depths (number of reads in each observed allele)" >',
#  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth" >',
#  '\n##FORMAT=<ID=FT,Number=1,Type=String,Description="Variant filters" >',
#  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Tumour ref count" >',
#  '\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Tumour alt counth" > \n'  
#      )
#      print(df)
# }
# prettyprint()



