

### analyze GISTIC output

## load required libraries
library("tidyverse")
library("plyranges")
## load GISTIC out put lesion file
all_lesions<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/all_lesions.conf_75.txt", header = TRUE)
all_lesions<-all_lesions[!duplicated(all_lesions$Descriptor),]

all_lesions$Wide.Peak.Limits<-gsub("\\s*\\([^\\)]+\\)","",all_lesions$Wide.Peak.Limits)

all_lesions_good<- all_lesions%>% separate(Wide.Peak.Limits, c("chr","start","end"))
all_lesions_good_selected<-all_lesions_good[,c("chr","start","end","Descriptor")]
all_lesions_good_selected$chr<-gsub("chr","",all_lesions_good_selected$chr)
gr_lesion <- makeGRangesFromDataFrame(all_lesions_good_selected, keep.extra.columns=TRUE)

### load cosmic cancer genes
cosmic_genehg19<-read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/cancer_genes/Census_May2023_hg19_fromCosmic.csv", sep = ",",header = TRUE)
cosmic_genehg19_good<-cosmic_genehg19%>%separate(Genome.Location,c("chrom","start","end"))


cosmic_genehg19_good<-cosmic_genehg19_good[,c("chrom","start","end","Gene.Symbol","Tier", "Hallmark")]

cosmic_genehg19_good<-cosmic_genehg19_good[cosmic_genehg19_good$start!="",]

gr_cosmic<-makeGRangesFromDataFrame(cosmic_genehg19_good, keep.extra.columns=TRUE)

findOverlaps(gr_lesion,gr_cosmic)
subsetByOverlaps(gr_lesion,gr_cosmic)


overlapped<-data.frame(join_overlap_inner(gr_lesion,gr_cosmic))
write.table(overlapped, file = "Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/overlapped_Gistic_with_Cosmochg1.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)