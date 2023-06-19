
### plot vcf stats

#https://krassowski.github.io/complex-upset/articles/Examples_R.html#fill-the-bars
###################################
##### analyze survivor consensus vcf files
###################################
library(ggplot2)
library(webr)
library(dplyr)
library(vcfR)
library(ggplot2)
library(tidyverse)
library(StructuralVariantAnnotation)
library("UpSetR")
library(ComplexUpset)
#library(stringi)


### => TO DO: drop contigs 
files_survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor", pattern = "*.vcf", full.names = TRUE,recursive=TRUE)

### lead info field
vcfs_survivor<-lapply(files_survivor,function(x){
     #readVcf(x)
     #info(readVcf(x))
     data.frame(info(readVcf(x)))
})


for (i in 1:length(vcfs_survivor)){
    vcfs_survivor[[i]]<-cbind(vcfs_survivor[[i]],files_survivor[i])
    }
type_data <- do.call("rbind", vcfs_survivor) 
rownames(type_data)<-NULL
names(type_data)[15]<-"path"

type_data$sample_id<-sub('.*/\\s*', '', type_data$path)
type_data$sample_id<-gsub("survivor_merged_filtered_","",gsub(".vcf","",type_data$sample_id))

type_data_good<-type_data%>%
        separate(sample_id,c("sample","technology","minimum_SVLEN","distance"))


type_data_good$sample <- gsub("1", "GCT",
        gsub("2", "RD",
        gsub("3", "SW982", type_data_good$sample)))


cols <- c("sample","technology","minimum_SVLEN","distance")
# create a new column `x` with the three columns collapsed together
type_data_good$uniq_ID <- apply( type_data_good[ , cols ] , 1 , paste , collapse = "_" )
#type_data_good<-type_data_good[,!names(type_data_good) %in% ("vcfs_combi.i.")]
     

type_data_good$SVLEN<-abs(type_data_good$SVLEN)
#data_plot<-type_data_good %>% 
#      mutate(SVCALLERS = sapply(SVCALLERS, toString))



type_data_good$length_range<-ifelse(type_data_good$SVLEN > 100 & type_data_good$SVLEN <= 300, "100-300",
                        ifelse(type_data_good$SVLEN > 300 & type_data_good$SVLEN <= 500, "300-500",
                        ifelse(type_data_good$SVLEN > 500 & type_data_good$SVLEN <= 1000, "500-1k",
                        ifelse(type_data_good$SVLEN > 1000 & type_data_good$SVLEN <= 10000, "1k-10k",
                        ifelse(type_data_good$SVLEN > 10000 & type_data_good$SVLEN <= 50000, "10k-50k",
                        ifelse(type_data_good$SVLEN > 50000, "> 50k",
                        "not_available"))))))





counts1 <- ddply(type_data_good, .(type_data_good$uniq_ID, type_data_good$SVTYPE, type_data_good$length_range), nrow)
colnames(counts1)<-c("sampleid_uniq","SVTYPE","length_range","freq")
write.table(counts1, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor/counts1_all_types.table", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

vcf_samples_count<-unique(counts1$sampleid_uniq)

for(ii in 1:length(vcf_samples_count)){

    focal_count<-counts1[counts1$sampleid_uniq==vcf_samples_count[ii],]
    #vec<-c("DEL","TRA","INS","DUP","INV")
    #counts1 <- counts1[order(factor(counts1$SVTYPE, levels=unique(vec))),]

    sample_names<-gsub("_","",vcf_samples_count[ii])
    pdf(file =paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor/piedonut_SVcount_",sample_names,".pdf"), height= 8, width = 6)
    PieDonut(focal_count, aes(SVTYPE , length_range, count=freq),r0=0.7,start=3*pi/2,labelpositionThreshold=0.1,showRatioThreshold = 0.000001,
         showRatioPie=TRUE,showRatioDonut=FALSE,labelposition=1,pieLabelSize=1.9,donutLabelSize=1.5,use.label=FALSE, title = paste(sample_names,"total number of variants",sum(focal_count$freq),sep = "_"))
     #print(plot_pie)
     dev.off()
    
        focal_count<-counts1[counts1$sampleid_uniq==vcf_samples_count[ii],]
        sample_names<-gsub("_","",vcf_samples_count[ii])

    pdf(file =paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor/applepie_SVcount_",sample_names,".pdf", sep = ""), height= 8, width = 6)
    plot_ptdo<-ggplot(focal_count, aes(x = SVTYPE, y = freq, fill = length_range)) +
    geom_col() +
    labs(title=paste(sample_names,"total number of variants",sum(focal_count$freq),sep = "_")) +
    scale_fill_viridis_d() +
    coord_polar("y")
    print(plot_ptdo)
    dev.off()

}



### barplot 
counts1_good<-counts1%>%
        separate(sampleid_uniq,c("sample","technology","minimum_SVLEN","distance"))

counts1_good$sample_distance<-paste(counts1_good$sample ,counts1_good$distance, sep = "_")
counts1_good$technology<-gsub("pacbio","PacBio",counts1_good$technology)
pdf(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor/barplot_count.pdf", height= 12, width=10)
ploty1<-ggplot(counts1_good,aes(x = technology,y=freq,fill=SVTYPE)) +
       geom_bar(stat="identity") +
       facet_grid(~sample_distance)+
       scale_fill_brewer(palette="Set1")
       #scale_fill_manual(values=c("red", "blue", "green")) #if you want to manually change 
print(ploty1)
dev.off()

############################################
############################################
##### amke upset :( plots of overlaps ######
############################################
############################################

### load matrix files 
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor", pattern = "*.txt", full.names = TRUE,recursive=TRUE)


matrix_tables<-lapply(matrix,function(x){
    read.csv(x, header = FALSE, strip.white=T, sep='')
})

for (i in 1:length(matrix_tables)){
    matrix_tables[[i]]<-cbind(matrix_tables[[i]],matrix[i])
    }
good_mat <- do.call("rbind", matrix_tables) 
names(good_mat)[5]<-"raw_id"


good_mat$sample_id<-sub('.*/\\s*', '', good_mat$raw_id)
good_mat$sample_id<-gsub("survivor_merged_filtered_","",gsub("_overlapped_3ways_matrix.txt","",good_mat$sample_id))

good_mat_fin<-good_mat%>%
        separate(sample_id,c("sample","technology","minimum_SVLEN","distance"))


good_mat_fin$sample <- gsub("1", "GCT",
        gsub("2", "RD",
        gsub("3", "SW982", type_data_good$sample)))

names(good_mat_fin)[1:4]<-c("cuteSV","nanoSV","sniffles","svim")
good_mat_fin<-good_mat_fin[,-5]



vcf_matrix<-cbind(type_data_good,good_mat_fin)
# columns to paste together
all_samples<-unique(vcf_matrix$uniq_ID)



for(jj in 1:length(all_samples)){

    vcf_matrix_focal<-vcf_matrix[vcf_matrix$uniq_ID ==all_samples[jj],]
    genres1 = colnames(vcf_matrix_focal)[22:25]
    vcf_matrix_focal[genres1] = vcf_matrix_focal[genres1] == 1
    t(head(vcf_matrix_focal[genres1], 3))
  

    pdf(file =paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor/upsetplot_SVcount_",all_samples[jj],".pdf", sep = ""), height= 12, width = 6)
    plot_upset<-upset(
    vcf_matrix_focal,
    genres1,
    base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='#ff4800'), guide='none')),
     mode='inclusive_intersection',
     annotations = list(
    # if not specified, the mode will follow the mode set in `upset()` call (here: `inclusive_intersection`)
    'Length (intersection size)'=(
      ggplot(mapping=aes(fill=SVTYPE))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
    ),
    'Length (SV length)'=(
      ggplot(mapping=aes(y=SVLEN))
      # checkout ggbeeswarm::geom_quasirandom for better results!
      + geom_jitter(aes(color=(SVLEN)), na.rm=TRUE)
      + geom_violin(alpha=0.5, na.rm=TRUE)

    )),
     min_size=1,
     width_ratio=0.1,n_intersections=6,
     set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))
  )
)
 #print(plot_upset)
 #dev.off()


### remove an outlier SV ####
    vcf_matrix_focal<-vcf_matrix[vcf_matrix$uniq_ID ==all_samples[jj],]
    genres1 = colnames(vcf_matrix_focal)[22:25]
    vcf_matrix_focal[genres1] = vcf_matrix_focal[genres1] == 1
    t(head(vcf_matrix_focal[genres1], 3))

    qq<-vcf_matrix_focal[vcf_matrix_focal$SVLEN!=max(vcf_matrix_focal$SVLEN),]
    qq<-qq[qq$SVLEN!=max(qq$SVLEN),]

    genresqq = colnames(qq)[22:25]
    qq[genresqq] = qq[genresqq] == 1
    t(head(qq[genresqq], 3))
  

    pdf(file =paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor/upsetplot_SVcount_log10SVLEN_drop2outlier",all_samples[jj],".pdf", sep = ""), height= 12, width = 6)
    plot_upset2<-upset(
    qq,
    genresqq,
    base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='#ff4800'), guide='none')),
     mode='inclusive_intersection',
     annotations = list(
    # if not specified, the mode will follow the mode set in `upset()` call (here: `inclusive_intersection`)
    'Length (intersection size)'=(
      ggplot(mapping=aes(fill=SVTYPE))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
    ),
    'Length (SV length)'=(
      ggplot(mapping=aes(y=SVLEN))
      # checkout ggbeeswarm::geom_quasirandom for better results!
      + geom_jitter(aes(color=log10(SVLEN)), na.rm=TRUE)
      + geom_violin(alpha=0.5, na.rm=TRUE)

    )),
     min_size=1,
     width_ratio=0.1,n_intersections=6,
     set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))
  )
)
print(plot_upset2)
 dev.off()

} ### main loop



### simplest upset plot
#set_size(8, 3)
#upset(vcf_matrix_focal, genres1, name='genre', width_ratio=0.1,n_intersections=6)  ## n_intersections=2

## with violin plot
#set_size(8, 8)
#set.seed(0)   # keep the same jitter for identical plots

## worked with violin plot
#upset(
#    vcf_matrix_focal,
#    genres1,
#    annotations = list(
#        
#        # 2nd method - using ggplot
#        'Rating'=(
#            # note that aes(x=intersection) is supplied by default and can be skipped
#            ggplot(mapping=aes(y=SVLEN))
#            # checkout ggbeeswarm::geom_quasirandom for better results!
#           + geom_jitter(aes(color=(SVLEN)), na.rm=TRUE)
#            + geom_violin(alpha=0.5, na.rm=TRUE)
#        )),
#    min_size=10,
#    width_ratio=0.1
#)

### worked stacked bar plot
#upset(
#    vcf_matrix_focal,
#    genres1,
#    annotations = list(
#        'MPAA Rating'=(
#            ggplot(mapping=aes(fill=SVTYPE))
#            + geom_bar(stat='count', position='fill')
#            + scale_y_continuous(labels=scales::percent_format())
#            
#            + ylab('MPAA Rating')
#        )
#    ),
#    width_ratio=0.1
#)



upset(
  vcf_matrix_focal,
  genres1,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='#ff4800'), guide='none')),
  mode='inclusive_intersection',
  annotations = list(
    # if not specified, the mode will follow the mode set in `upset()` call (here: `inclusive_intersection`)
    'Length (intersection size)'=(
      ggplot(mapping=aes(fill=SVTYPE))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
    ),
    'Length (SV length)'=(
      ggplot(mapping=aes(y=SVLEN))
      # checkout ggbeeswarm::geom_quasirandom for better results!
      + geom_jitter(aes(color=(SVLEN)), na.rm=TRUE)
      + geom_violin(alpha=0.5, na.rm=TRUE)

    )),
  min_size=1,
  width_ratio=0.1,n_intersections=6,
  set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))
)
)























vcfs_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio1/sv_calling/survivor", pattern = "*.vcf", full.names = TRUE,recursive=TRUE)
### lead info field
vcfs_consensus<-lapply(vcfs_file,function(x){
     #readVcf(x)
     #info(readVcf(x))
     data.frame(info(readVcf(x)))
})


for (i in 1:length(vcfs_consensus)){
    vcfs_consensus[[i]]<-cbind(vcfs_consensus[[i]],vcfs_file[i])
    }

vcfs_consensus_good<- do.call("rbind", vcfs_consensus) 
names(vcfs_consensus_good)[15]<-"raw_id"

vcfs_consensus_good$sample_id<-sub('.*/\\s*', '', vcfs_consensus_good$raw_id)
vcfs_consensus_good$sample_id<-gsub("survivor_merged_filtered_","",gsub(".vcf","",vcfs_consensus_good$sample_id))

final_consensus<-vcfs_consensus_good%>%
        separate(sample_id,c("sample","technology","minimum_SVLEN","distance"))


final_consensus$sample <- gsub("1", "GCT",
        gsub("2", "RD",
        gsub("3", "SW982", final_consensus$sample)))

final_consensus<-final_consensus[,!names(final_consensus) %in% ("vcfs_file[i]")]
      


vcf_matrix<-cbind(final_consensus,good_mat_fin)
# columns to paste together
cols <- c("sample","technology","minimum_SVLEN","distance")

# create a new column `x` with the three columns collapsed together
vcf_matrix$uniq_ID <- apply( vcf_matrix[ , cols ] , 1 , paste , collapse = "_" )
vcf_matrix$SVLEN<-abs(vcf_matrix$SVLEN)



vcf_matrix_focal<-vcf_matrix[vcf_matrix$uniq_ID =="GCT_pacbio_SVLEN100_DIS500",]

vcf_matrix_focal[genres1] = vcf_matrix_focal[genres1] == 1
t(head(vcf_matrix_focal[genres1], 3))

### simplest upset plot
#set_size(8, 3)
#upset(vcf_matrix_focal, genres1, name='genre', width_ratio=0.1,n_intersections=6)  ## n_intersections=2

## with violin plot
#set_size(8, 8)
#set.seed(0)   # keep the same jitter for identical plots

## worked with violin plot
#upset(
#    vcf_matrix_focal,
#    genres1,
#    annotations = list(
#        
#        # 2nd method - using ggplot
#        'Rating'=(
#            # note that aes(x=intersection) is supplied by default and can be skipped
#            ggplot(mapping=aes(y=SVLEN))
#            # checkout ggbeeswarm::geom_quasirandom for better results!
#           + geom_jitter(aes(color=(SVLEN)), na.rm=TRUE)
#            + geom_violin(alpha=0.5, na.rm=TRUE)
#        )),
#    min_size=10,
#    width_ratio=0.1
#)

### worked stacked bar plot
#upset(
#    vcf_matrix_focal,
#    genres1,
#    annotations = list(
#        'MPAA Rating'=(
#            ggplot(mapping=aes(fill=SVTYPE))
#            + geom_bar(stat='count', position='fill')
#            + scale_y_continuous(labels=scales::percent_format())
#            
#            + ylab('MPAA Rating')
#        )
#    ),
#    width_ratio=0.1
#)



upset(
  vcf_matrix_focal,
  genres1,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill='bars_color')
    ) + scale_fill_manual(values=c('bars_color'='#ff4800'), guide='none')),
  mode='inclusive_intersection',
  annotations = list(
    # if not specified, the mode will follow the mode set in `upset()` call (here: `inclusive_intersection`)
    'Length (intersection size)'=(
      ggplot(mapping=aes(fill=SVTYPE))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
    ),
    'Length (SV length)'=(
      ggplot(mapping=aes(y=SVLEN))
      # checkout ggbeeswarm::geom_quasirandom for better results!
      + geom_jitter(aes(color=(SVLEN)), na.rm=TRUE)
      + geom_violin(alpha=0.5, na.rm=TRUE)

    )),
  min_size=1,
  width_ratio=0.1,n_intersections=6,
  set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))
)
)


#########################
#### venn diagram in r
#########################
library(VennDiagram)
#t=read.table("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/nanopore1/sv_calling/survivor/survivor_3matched_SVLEN100_DIS500/3_merged_overlapp.txt",header=F)

cols <- c("sample","technology","minimum_SVLEN","distance")

# create a new column `x` with the three columns collapsed together
good_mat_fin$uniq_ID <- apply( good_mat_fin[ , cols ] , 1 , paste , collapse = "_" )


venn_samples<-unique(good_mat_fin$uniq_ID)
venn_final<-good_mat_fin[good_mat_fin$uniq_ID=="GCT_pacbio_SVLEN100_DIS500",]

venn.diagram(list
             (cuteSV=which(venn_final[,1]==1), nanoSV=which(venn_final[,2]==1), Sniffles=which(venn_final[,3]==1),svim=which(venn_final[,4]==1)) , 
             fill = c("green", "orange", "blue","tomato") , 
             alpha = c(0.5, 0.5, 0.5, 0.5), 
             cex = 2, lty =2, ext.text=TRUE, 
             filename = "~/Desktop/3_ONT.jpeg")

