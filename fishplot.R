
#### make fish tree plots for pyclone outputs
rm(list = ls())
#library(fishplot)
library(tidyverse)
library(clonevol)
library("plyr")

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/cfDNA_analysis/Updated_files_pyclone/")


filenames_py <- list.files("deepsequencing_purity_full",pattern="*.pyclone.results.tsv", full.names = TRUE)
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

py_data <- do.call("rbind", attackStats_py) 
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))
#py_data_good<-py_data_good[!grepl("_C",py_data_good$sample_id),]


py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")


meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)
#paied_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167")
#meta_paied<-meta[meta$unique_sample_id%in%paied_samples,]

py_samples<-unique(py_data_good$unique_sample_id)
#[1] "SRC125"  "SRC127"  "SRC130"  "SRC150"  "SRC167"  "SRC168"  "SRC169"  "SRC170"  "SRC171" 
#[10] "SRC172"  "SRC173"  "TB11985" "TB12052" "TB13092" "TB13712" "TB13959" "TB22446" "TB8016" 
#[19] "TB9051"  "TB9573" 

################
#out_res_py<-NULL
#for (zz in 1:length(py_samples)){

#py_focal<-py_data_good[py_data_good$unique_sample_id==py_samples[zz],]
py_focal<-py_data_good[py_data_good$unique_sample_id=="SRC169",]
mutaion_RT_focal<-merge(py_focal,meta, by.x = "sample_id", by.y = "sample_id")
mutaion_RT_focal$identifier<-paste(mutaion_RT_focal$unique_sample_id.x,mutaion_RT_focal$RT_status, sep = "_")
mutaion_RT_focal$code<-paste(mutaion_RT_focal$sample_id,mutaion_RT_focal$RT_status, sep = "_")
mutaion_RT_focal$pos_identifier<-gsub("chr","",mutaion_RT_focal$pos_id)

pos_identifier<-unique(mutaion_RT_focal$pos_identifier)
mutaion_RT_focal_good<-mutaion_RT_focal[,c("code","cluster_id","cellular_prevalence","pos_identifier")]
mutaion_RT_focal_good$code<-gsub("-.*","",mutaion_RT_focal_good$code)
#mutaion_RT_focal_good<-mutaion_RT_focal_good[!(grepl("SRC168_7_noRT",mutaion_RT_focal_good$code) |  grepl("SRC168_8_noRT",mutaion_RT_focal_good$code) | grepl("SRC168_9_noRT",mutaion_RT_focal_good$cod)),]


out_res_phy<-NULL
for(ii in 1:length(pos_identifier)) {

    focal_pos<-mutaion_RT_focal_good[mutaion_RT_focal_good$pos_identifier== pos_identifier[ii],]
    focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)*100
    #focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)/2

    aa<-data.frame(t(focal_pos))
    colnames(aa)<-aa[1,]
    bb<-aa[3,]
    rownames(bb)<-NULL
    colnames(bb)<-paste(colnames(bb),"cp", sep = "_")
    bb[] <- lapply(bb, function(x) as.numeric(as.character(x)))
        
    bb<-bb/2
    colnames(bb)<-gsub("cp","vaf",colnames(bb))
    #both<-cbind(bb,vaf)
    bb$cluster<-unique(focal_pos$cluster_id)
    bb$position<-unique(focal_pos$pos_identifier)
    #both$sample<-py_samples[zz]
    #both$sample<-"TB9051"
    out_res_phy<-rbind(bb,out_res_phy)

          #data_clon<-data.frame(t(aggregate( cellular_prevalence ~ RT_status, focal_pos, mean )))
          #names(data_clon)<-NULL
          #data_clon<-data_clon[-1,]

          #rownames(data_clon)<-NULL
          #colnames(data_clon)<-c("postRT.cp","pretRT.cp")
          #data_clon$cluster<-unique(focal_pos$cluster_id)
          #data_clon$pos_id<-unique(focal_pos$pos_identifier)
          #out_res_phy<-rbind(data_clon,out_res_phy)
       }
#}



#### take care of founding population
vaf.col.names = grep('.vaf', colnames(out_res_phy), value=T)
colnames(out_res_phy) = gsub('.vaf', '', colnames(out_res_phy))
vaf.col.names = gsub('.vaf', '', vaf.col.names)


#ccf.col.names = grep('.cp', colnames(out_res_phy), value=T)
#colnames(out_res_phy) = gsub('.cp', '', colnames(out_res_phy))
#ccf.col.names = gsub('.cp', '', vaf.col.names)

# make sure cluster are continuous integer starting at 1
out_res_phy$cluster[out_res_phy$cluster == 0] = max(out_res_phy$cluster) + 1


out_res_phy$cluster[out_res_phy$cluster == 4] <- 77
out_res_phy$cluster[out_res_phy$cluster == 1] <- 4
out_res_phy$cluster[out_res_phy$cluster == 77] <- 1



# let clonevol decide colors
clone.colors <-NULL

# plot cluster of variants to make sure clustering is good
pdf(file = paste("~/Desktop/clonevol.vaf.box_TB9051.pdf", sep = "") ,width = 5, height = 10, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(out_res_phy,
    cluster.col.name = 'cluster',
    show.cluster.size = FALSE,
    cluster.size.text.color = 'tomato',
    vaf.col.names = vaf.col.names,
    vaf.limits = 70,
    sample.title.size = 4,
    violin = FALSE,
    box = FALSE,
    jitter = TRUE,
    jitter.shape = 1,
    jitter.color = clone.colors,
    jitter.size = 0.8,
    jitter.alpha = 1,
    jitter.center.method = 'median',
    jitter.center.size = 1,
    jitter.center.color = 'darkgray',
    jitter.center.display.value = 'none',
    highlight = 'is.driver',
    highlight.shape = 21,
    highlight.color = 'tomato',
    highlight.fill.color = 'green',
    highlight.note.col.name = 'gene',
    highlight.note.size = 1,
    order.by.total.vaf = FALSE)
    #print(pp)
    dev.off()

### plot paiwise VAF 

clone.colors = NULL
#pdf(file = paste("deepsequencing_purity_full/clonevol_outs/clonevol.pairwise.vaf.box_",py_samples[zz],".pdf", sep = "") ,width = 5, height = 5, useDingbats = FALSE, title='')
#plot.pairwise(out_res_phy, col.names = vaf.col.names,
#out.prefix = 'variants.pairwise.plot',
#colors = clone.colors)
#dev.off()

clone.colors = NULL
# Plotting mean/median of clusters across samples (cluster flow)
#pdf(file = paste("deepsequencing_purity_full/clonevol_outs/clonevol.vaf.across.samples.box_",py_samples[zz],".pdf", sep = "") ,width = 15, height = 5, useDingbats = FALSE, title='')
pdf(file = "~/Desktop/clonevol.vaf.across.samples.box_SRC169.pdf", width = 15, height = 5, useDingbats = FALSE)
plot.cluster.flow(out_res_phy, vaf.col.names = vaf.col.names,
sample.names = vaf.col.names,
colors = clone.colors)
dev.off()


y = infer.clonal.models(variants = out_res_phy,
    cluster.col.name = 'cluster',
    vaf.col.names = vaf.col.names,
    sample.groups = NULL,
    cancer.initiation.model='monoclonal',
    subclonal.test = 'bootstrap',
    subclonal.test.model = 'non-parametric',
    num.boots = 1000,
    founding.cluster = 1,
    cluster.center = 'mean',
    ignore.clusters = NULL,
    clone.colors = NULL,
    min.cluster.vaf = 0.005,
    sum.p = 0.05,
    alpha = 0.05)
    
y = convert.consensus.tree.clone.to.branch(y)   



# now plot
plot.clonal.models(y,
    box.plot = TRUE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.highlight = 'is.driver',
    fancy.variant.boxplot.highlight.shape = 21,
    fancy.variant.boxplot.highlight.fill.color = 'red',
    fancy.variant.boxplot.highlight.color = 'black',
    fancy.variant.boxplot.highlight.note.col.name = 'gene',
    fancy.variant.boxplot.highlight.note.color = 'blue',
    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1,
    fancy.variant.boxplot.jitter.center.color = 'grey50',
    fancy.variant.boxplot.base_size = 12,
    fancy.variant.boxplot.plot.margin = 1,
    fancy.variant.boxplot.vaf.suffix = '.VAF',
    clone.shape = 'bell',
    bell.event = TRUE,
    bell.event.label.color = 'blue',
    bell.event.label.angle = 60,
    clone.time.step.scale = 1,
    bell.curve.step = 2,
    merged.tree.plot = TRUE,
    tree.node.label.split.character = NULL,
    tree.node.shape = 'circle',
    tree.node.size = 30,
    tree.node.text.size = 0.5,
    merged.tree.node.size.scale = 1.25,
    merged.tree.node.text.size.scale = 2.5,
    merged.tree.cell.frac.ci = FALSE,
    merged.tree.clone.as.branch = TRUE,
    mtcab.event.sep.char = ',',
    mtcab.branch.text.size = 1,
    mtcab.branch.width = 1,
    mtcab.node.size = 3,
    mtcab.node.label.size = 1,
    mtcab.node.text.size = 1.5,
    cell.plot = TRUE, num.cells = 100,
    cell.border.size = 0.25,
    cell.border.color = 'black',
    clone.grouping = 'horizontal',
    scale.monoclonal.cell.frac = TRUE,
    show.score = FALSE,
    cell.frac.ci = TRUE,
    disable.cell.frac = FALSE,
    out.dir = 'output_TB13092_regional', out.format = 'pdf',
    overwrite.output = TRUE,
    width = 11, height = 7,
    panel.widths = c(3,4,2,4,4))


f = generateFishplotInputs(results=y)
fishes = createFishPlotObjects(f)
#plot with fishplot
pdf(file = "fishplot.src125.notmerged.reagional.pdf", width=6, height=5)
for (i in 1:1){
fish = layoutClones(fishes[[i]])
fish = setCol(fish,f$clonevol.clone.colors)
fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
vlines=seq(1, length(vaf.col.names)), vlab=vaf.col.names, pad.left=0.5)
}
dev.off()



#####################################################
####### calulate mean for preRT and postRT samples 


limit_to_vaf<-out_res_phy[,vaf.col.names]
preRT<-out_res_phy[,grep("pre",names(limit_to_vaf))]
#preRT<-data.frame("SRC167_5_preRT"=(limit_to_vaf[,grep("pre",names(limit_to_vaf))]))
postRT<-limit_to_vaf[,grep("post",names(limit_to_vaf))]

#preRT<-out_res_phy[,c(1:6)]
#postRT<-out_res_phy[,c(7:8)]

preRT_mean<-data.frame("pre.vaf"=rowMeans(preRT))
postRT_mean<-data.frame("post.vaf"=rowMeans(postRT))
cluste<-out_res_phy[,c("cluster","position")]
both<-cbind(preRT_mean,postRT_mean,cluste)



vaf.col.names = grep('.vaf', colnames(both), value=T)
colnames(both) = gsub('.vaf', '', colnames(both))
vaf.col.names = gsub('.vaf', '', vaf.col.names)


#both$cluster[both$cluster == 2] <- 77
#both$cluster[both$cluster == 1] <- 2
#both$cluster[both$cluster == 77] <- 1

# make sure cluster are continuous integer starting at 1
#both$cluster[both$cluster == 0] = max(both$cluster) + 1

# let clonevol decide colors
clone.colors <-NULL

### plot variant cluster
pdf(file = "~/Desktop/clonevol.vaf_mean_regional.box_SRC169.pdf" ,width = 5, height = 4, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(both,
    cluster.col.name = 'cluster',
    show.cluster.size = FALSE,
    cluster.size.text.color = 'blue',
    vaf.col.names = vaf.col.names,
    vaf.limits = 100,
    sample.title.size = 20,
    violin = FALSE,
    box = FALSE,
    jitter = TRUE,
    jitter.shape = 1,
    jitter.color = clone.colors,
    jitter.size = 1,
    jitter.alpha = 1,
    jitter.center.method = 'median',
    jitter.center.size = 1,
    jitter.center.color = 'darkgray',
    jitter.center.display.value = 'none',
    highlight = 'is.driver',
    highlight.shape = 21,
    highlight.color = 'blue',
    highlight.fill.color = 'green',
    highlight.note.col.name = 'gene',
    highlight.note.size = 2,
    order.by.total.vaf = FALSE)
    dev.off()

### plot paiwise VAF 
#plot.pairwise(both, col.names = vaf.col.names,
#out.prefix = 'variants.pairwise.plot',
#colors = clone.colors)

clone.colors<-NULL
# Plotting mean/median of clusters across samples (cluster flow)
pdf(file = "~/Desktop/clonevol.vaf_mean_regional.across.samples.box_SRC127.pdf",width = 5, height = 5, useDingbats = FALSE, title='')
plot.cluster.flow(both, vaf.col.names = vaf.col.names,
sample.names = c("preRT","postRT"),
colors = clone.colors)
dev.off()

y = infer.clonal.models(variants = both,
    cluster.col.name = 'cluster',
    vaf.col.names = vaf.col.names,
    sample.groups = NULL,
    cancer.initiation.model='monoclonal',
    subclonal.test = 'bootstrap',
    subclonal.test.model = 'non-parametric',
    num.boots = 1000,
    founding.cluster = 1,
    cluster.center = 'mean',
    ignore.clusters = NULL,
    clone.colors = NULL,
    min.cluster.vaf = 0.00,
    sum.p = 0.05,
    alpha = 0.05)


y = convert.consensus.tree.clone.to.branch(y)

# now plot
clone.colors<-NULL
plot.clonal.models(y,
    box.plot = TRUE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.highlight = 'is.driver',
    fancy.variant.boxplot.highlight.shape = 21,
    fancy.variant.boxplot.highlight.fill.color = 'red',
    fancy.variant.boxplot.highlight.color = 'black',
    fancy.variant.boxplot.highlight.note.col.name = 'gene',
    fancy.variant.boxplot.highlight.note.color = 'blue',
    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1,
    fancy.variant.boxplot.jitter.center.color = 'grey50',
    fancy.variant.boxplot.base_size = 12,
    fancy.variant.boxplot.plot.margin = 1,
    fancy.variant.boxplot.vaf.suffix = '.VAF',
    clone.shape = 'bell',
    bell.event = TRUE,
    bell.event.label.color = 'blue',
    bell.event.label.angle = 60,
    clone.time.step.scale = 1,
    bell.curve.step = 2,
    merged.tree.plot = TRUE,
    tree.node.label.split.character = NULL,
    tree.node.shape = 'circle',
    tree.node.size = 30,
    tree.node.text.size = 0.5,
    merged.tree.node.size.scale = 1.25,
    merged.tree.node.text.size.scale = 2.5,
    merged.tree.cell.frac.ci = FALSE,
    merged.tree.clone.as.branch = TRUE,
    mtcab.event.sep.char = ',',
    mtcab.branch.text.size = 1,
    mtcab.branch.width = 0.3,
    mtcab.node.size = 3,
    mtcab.node.label.size = 1,
    mtcab.node.text.size = 1.5,
    cell.plot = TRUE, num.cells = 100,
    cell.border.size = 0.25,
    cell.border.color = 'black',
    clone.grouping = 'horizontal',
    scale.monoclonal.cell.frac = TRUE,
    show.score = FALSE,
    cell.frac.ci = TRUE,
    disable.cell.frac = FALSE,
    out.dir = '~/Desktop/output.src127', out.format = 'pdf',
    overwrite.output = TRUE,
    width = 11, height = 7,
    panel.widths = c(3,4,2,4,4))

### fishy
f = generateFishplotInputs(results=y)
fishes = createFishPlotObjects(f)
#plot with fishplot
pdf(file = "~/Desktop/fishplot.src127.mean.reagional.pdf", width=6, height=5)
for (i in 1:2){
fish = layoutClones(fishes[[i]])
fish = setCol(fish,f$clonevol.clone.colors)
fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
vlines=seq(1, length(vaf.col.names)), vlab=vaf.col.names, pad.left=0.5)
}
dev.off()




############################################
############################################
######## Single Time Point Samples #########
############################################
############################################

rm(list = ls())
#library(fishplot)
library(tidyverse)
library(clonevol)
library("plyr")

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/cfDNA_analysis/Updated_files_pyclone/")


filenames_py <- list.files("deepsequencing_purity_full",pattern="*.pyclone.results.tsv", full.names = TRUE)
attackStats_py <- lapply(filenames_py,function(x) {
     read.delim(x, header=TRUE, sep = "\t")
     })

py_data <- do.call("rbind", attackStats_py) 
py_data_good<-py_data %>% separate(mutation_id, c('Chrom', 'Pos','Ref' ,'Alt' ))
#py_data_good<-py_data_good[!grepl("_C",py_data_good$sample_id),]


py_data_good$unique_sample_id<-gsub("_.*$","",py_data_good$sample_id)
py_data_good$pos_id<-paste(py_data_good$Chrom,py_data_good$Pos, sep = "_")


meta<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_inputs/metadata_updated.txt", header = FALSE, sep= "\t")[,c(1:4,6)]
colnames(meta)<-c("sample_id", "purity","ploidy","RT_status","RT_code")
meta$unique_sample_id<-gsub("_.*$","",meta$sample_id)
#paied_samples<-c("SRC125","SRC127","SRC169","SRC170","SRC171","TB9051","SRC167")
#meta_paied<-meta[meta$unique_sample_id%in%paied_samples,]

py_samples<-unique(py_data_good$unique_sample_id)
#[1] "SRC125"  "SRC127"  "SRC130"  "SRC150"  "SRC167"  "SRC168"  "SRC169"  "SRC170"  "SRC171" 
#[10] "SRC172"  "SRC173"  "TB11985" "TB12052" "TB13092" "TB13712" "TB13959" "TB22446" "TB8016" 
#[19] "TB9051"  "TB9573" 

################
#out_res_py<-NULL
#for (zz in 1:length(py_samples)){
#py_focal<-py_data_good[py_data_good$unique_sample_id==py_samples[zz],]
py_focal<-py_data_good[py_data_good$unique_sample_id=="SRC130",]
mutaion_RT_focal<-merge(py_focal,meta, by.x = "sample_id", by.y = "sample_id")
mutaion_RT_focal$identifier<-paste(mutaion_RT_focal$unique_sample_id.x,mutaion_RT_focal$RT_status, sep = "_")
mutaion_RT_focal$code<-paste(mutaion_RT_focal$sample_id,mutaion_RT_focal$RT_status, sep = "_")
mutaion_RT_focal$pos_identifier<-gsub("chr","",mutaion_RT_focal$pos_id)

pos_identifier<-unique(mutaion_RT_focal$pos_identifier)
mutaion_RT_focal_good<-mutaion_RT_focal[,c("code","cluster_id","cellular_prevalence","pos_identifier")]
mutaion_RT_focal_good$code<-gsub("-.*","",mutaion_RT_focal_good$code)
#mutaion_RT_focal_good<-mutaion_RT_focal_good[!(grepl("SRC172_6_noRT",mutaion_RT_focal_good$code) |  grepl("SRC172_5_noRT",mutaion_RT_focal_good$code)),]


out_res_phy<-NULL
for(ii in 1:length(pos_identifier)) {

    focal_pos<-mutaion_RT_focal_good[mutaion_RT_focal_good$pos_identifier== pos_identifier[ii],]
    focal_pos$cellular_prevalence<-(focal_pos$cellular_prevalence)*100
    aa<-data.frame(t(focal_pos))
    colnames(aa)<-aa[1,]
    bb<-aa[3,]
    rownames(bb)<-NULL
    colnames(bb)<-paste(colnames(bb),"cp", sep = "_")
    bb[] <- lapply(bb, function(x) as.numeric(as.character(x)))
        
    vaf<-bb/2
    colnames(vaf)<-gsub("cp","vaf",colnames(vaf))
    both<-cbind(bb,vaf)
    both$cluster<-unique(focal_pos$cluster_id)
    both$position<-unique(focal_pos$pos_identifier)
    #both$sample<-py_samples[zz]
    #both$sample<-"TB9051"
    out_res_phy<-rbind(both,out_res_phy)

          #data_clon<-data.frame(t(aggregate( cellular_prevalence ~ RT_status, focal_pos, mean )))
          #names(data_clon)<-NULL
          #data_clon<-data_clon[-1,]

          #rownames(data_clon)<-NULL
          #colnames(data_clon)<-c("postRT.cp","pretRT.cp")
          #data_clon$cluster<-unique(focal_pos$cluster_id)
          #data_clon$pos_id<-unique(focal_pos$pos_identifier)
          #out_res_phy<-rbind(data_clon,out_res_phy)
       }
#}


#### take care of founding population
vaf.col.names = grep('.vaf', colnames(out_res_phy), value=T)
colnames(out_res_phy) = gsub('.vaf', '', colnames(out_res_phy))
vaf.col.names = gsub('.vaf', '', vaf.col.names)


#ccf.col.names = grep('.cp', colnames(out_res_phy), value=T)
#colnames(out_res_phy) = gsub('.cp', '', colnames(out_res_phy))
#ccf.col.names = gsub('.cp', '', vaf.col.names)

# make sure cluster are continuous integer starting at 1
out_res_phy$cluster[out_res_phy$cluster == 0] = max(out_res_phy$cluster) + 1


out_res_phy$cluster[out_res_phy$cluster == 2] <- 77
out_res_phy$cluster[out_res_phy$cluster == 1] <- 2
out_res_phy$cluster[out_res_phy$cluster == 77] <- 1



# let clonevol decide colors
clone.colors <-NULL

# plot cluster of variants to make sure clustering is good
pdf(file = paste("~/Desktop/clonevol.vaf.box_TB22446.pdf", sep = "") ,width = 5, height = 10, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(out_res_phy,
    cluster.col.name = 'cluster',
    show.cluster.size = FALSE,
    cluster.size.text.color = 'tomato',
    vaf.col.names = vaf.col.names,
    vaf.limits = 70,
    sample.title.size = 4,
    violin = FALSE,
    box = FALSE,
    jitter = TRUE,
    jitter.shape = 1,
    jitter.color = clone.colors,
    jitter.size = 0.8,
    jitter.alpha = 1,
    jitter.center.method = 'median',
    jitter.center.size = 1,
    jitter.center.color = 'darkgray',
    jitter.center.display.value = 'none',
    highlight = 'is.driver',
    highlight.shape = 21,
    highlight.color = 'tomato',
    highlight.fill.color = 'green',
    highlight.note.col.name = 'gene',
    highlight.note.size = 1,
    order.by.total.vaf = FALSE)
    #print(pp)
    dev.off()

### plot paiwise VAF 

clone.colors = NULL
#pdf(file = paste("deepsequencing_purity_full/clonevol_outs/clonevol.pairwise.vaf.box_",py_samples[zz],".pdf", sep = "") ,width = 5, height = 5, useDingbats = FALSE, title='')
#plot.pairwise(out_res_phy, col.names = vaf.col.names,
#out.prefix = 'variants.pairwise.plot',
#colors = clone.colors)
#dev.off()

clone.colors = NULL
# Plotting mean/median of clusters across samples (cluster flow)
#pdf(file = paste("deepsequencing_purity_full/clonevol_outs/clonevol.vaf.across.samples.box_",py_samples[zz],".pdf", sep = "") ,width = 15, height = 5, useDingbats = FALSE, title='')
pdf(file = "~/Desktop/clonevol.vaf.across.samples.box_SRC130.pdf", width = 15, height = 5, useDingbats = FALSE)
plot.cluster.flow(out_res_phy, vaf.col.names = vaf.col.names,
sample.names = vaf.col.names,
colors = clone.colors)
dev.off()


y = infer.clonal.models(variants = out_res_phy,
    cluster.col.name = 'cluster',
    vaf.col.names = vaf.col.names,
    sample.groups = NULL,
    cancer.initiation.model='monoclonal',
    subclonal.test = 'bootstrap',
    subclonal.test.model = 'non-parametric',
    num.boots = 1000,
    founding.cluster = 1,
    cluster.center = 'mean',
    ignore.clusters = NULL,
    clone.colors = NULL,
    min.cluster.vaf = 0.00,
    sum.p = 0.05,
    alpha = 0.05)
    
y = convert.consensus.tree.clone.to.branch(y)   



# now plot
plot.clonal.models(y,
    box.plot = TRUE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.highlight = 'is.driver',
    fancy.variant.boxplot.highlight.shape = 21,
    fancy.variant.boxplot.highlight.fill.color = 'red',
    fancy.variant.boxplot.highlight.color = 'black',
    fancy.variant.boxplot.highlight.note.col.name = 'gene',
    fancy.variant.boxplot.highlight.note.color = 'blue',
    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1,
    fancy.variant.boxplot.jitter.center.color = 'grey50',
    fancy.variant.boxplot.base_size = 12,
    fancy.variant.boxplot.plot.margin = 1,
    fancy.variant.boxplot.vaf.suffix = '.VAF',
    clone.shape = 'bell',
    bell.event = TRUE,
    bell.event.label.color = 'blue',
    bell.event.label.angle = 60,
    clone.time.step.scale = 1,
    bell.curve.step = 2,
    merged.tree.plot = TRUE,
    tree.node.label.split.character = NULL,
    tree.node.shape = 'circle',
    tree.node.size = 30,
    tree.node.text.size = 0.5,
    merged.tree.node.size.scale = 1.25,
    merged.tree.node.text.size.scale = 2.5,
    merged.tree.cell.frac.ci = FALSE,
    merged.tree.clone.as.branch = TRUE,
    mtcab.event.sep.char = ',',
    mtcab.branch.text.size = 1,
    mtcab.branch.width = 1,
    mtcab.node.size = 3,
    mtcab.node.label.size = 1,
    mtcab.node.text.size = 1.5,
    cell.plot = TRUE, num.cells = 100,
    cell.border.size = 0.25,
    cell.border.color = 'black',
    clone.grouping = 'horizontal',
    scale.monoclonal.cell.frac = TRUE,
    show.score = FALSE,
    cell.frac.ci = TRUE,
    disable.cell.frac = FALSE,
    out.dir = 'output_SRC168_notmerged_3samples_removed', out.format = 'pdf',
    overwrite.output = TRUE,
    width = 11, height = 7,
    panel.widths = c(3,4,2,4,4))


f = generateFishplotInputs(results=y)
fishes = createFishPlotObjects(f)
#plot with fishplot
pdf(file = "fishplot.src125.notmerged.reagional.pdf", width=6, height=5)
for (i in 1:1){
fish = layoutClones(fishes[[i]])
fish = setCol(fish,f$clonevol.clone.colors)
fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
vlines=seq(1, length(vaf.col.names)), vlab=vaf.col.names, pad.left=0.5)
}
dev.off()


#############################################################
#############################################################
####### calulate mean for preRT and postRT samples #########

#### take care of founding population
vaf.col.names = grep('.vaf', colnames(out_res_phy), value=T)
colnames(out_res_phy) = gsub('.vaf', '', colnames(out_res_phy))
vaf.col.names = gsub('.vaf', '', vaf.col.names)


#ccf.col.names = grep('.cp', colnames(out_res_phy), value=T)
#colnames(out_res_phy) = gsub('.cp', '', colnames(out_res_phy))
#ccf.col.names = gsub('.cp', '', vaf.col.names)

# make sure cluster are continuous integer starting at 1
out_res_phy$cluster[out_res_phy$cluster == 0] = max(out_res_phy$cluster) + 1


out_res_phy$cluster[out_res_phy$cluster == 6] <- 77
out_res_phy$cluster[out_res_phy$cluster == 1] <- 6
out_res_phy$cluster[out_res_phy$cluster == 77] <- 1



limit_to_vaf<-out_res_phy[,vaf.col.names]
noRT_mean<-data.frame("noRT.vaf"=rowMeans(limit_to_vaf))

dummy<-noRT_mean
colnames(dummy)<-"dummy.vaf"

cluste<-out_res_phy[,c("cluster","position")]
both<-cbind(noRT_mean,dummy,cluste)
vaf.col.names<-c("noRT.vaf","dummy.vaf")

# make sure cluster are continuous integer starting at 1
#both$cluster[both$cluster == 0] = max(both$cluster) + 1


#both$cluster[both$cluster == 2] <- 77
#both$cluster[both$cluster == 1] <- 2
#both$cluster[both$cluster == 77] <- 1

# make sure cluster are continuous integer starting at 1
#both$cluster[both$cluster == 0] = max(both$cluster) + 1

# let clonevol decide colors
clone.colors <-NULL

### plot variant cluster
pdf(file = "~/Desktop/clonevol.vaf_mean_regional.box_TB13092.pdf" ,width = 5, height = 4, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(both,
    cluster.col.name = 'cluster',
    show.cluster.size = FALSE,
    cluster.size.text.color = 'blue',
    vaf.col.names = vaf.col.names,
    vaf.limits = 100,
    sample.title.size = 20,
    violin = FALSE,
    box = FALSE,
    jitter = TRUE,
    jitter.shape = 1,
    jitter.color = clone.colors,
    jitter.size = 1,
    jitter.alpha = 1,
    jitter.center.method = 'median',
    jitter.center.size = 1,
    jitter.center.color = 'darkgray',
    jitter.center.display.value = 'none',
    highlight = 'is.driver',
    highlight.shape = 21,
    highlight.color = 'blue',
    highlight.fill.color = 'green',
    highlight.note.col.name = 'gene',
    highlight.note.size = 2,
    order.by.total.vaf = FALSE)
    dev.off()

y = infer.clonal.models(variants = both,
    cluster.col.name = 'cluster',
    vaf.col.names = vaf.col.names,
    sample.groups = NULL,
    cancer.initiation.model='monoclonal',
    subclonal.test = 'bootstrap',
    subclonal.test.model = 'non-parametric',
    num.boots = 1000,
    founding.cluster = 1,
    cluster.center = 'mean',
    ignore.clusters = NULL,
    clone.colors = NULL,
    min.cluster.vaf = 0.000,
    sum.p = 0.05,
    alpha = 0.05)


y = convert.consensus.tree.clone.to.branch(y)

# now plot
clone.colors<-NULL
plot.clonal.models(y,
    box.plot = TRUE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.highlight = 'is.driver',
    fancy.variant.boxplot.highlight.shape = 21,
    fancy.variant.boxplot.highlight.fill.color = 'red',
    fancy.variant.boxplot.highlight.color = 'black',
    fancy.variant.boxplot.highlight.note.col.name = 'gene',
    fancy.variant.boxplot.highlight.note.color = 'blue',
    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1,
    fancy.variant.boxplot.jitter.center.color = 'grey50',
    fancy.variant.boxplot.base_size = 12,
    fancy.variant.boxplot.plot.margin = 1,
    fancy.variant.boxplot.vaf.suffix = '.VAF',
    clone.shape = 'bell',
    bell.event = TRUE,
    #bell.event.label.color = 'blue',
    bell.event.label.angle = 60,
    clone.time.step.scale = 1,
    bell.curve.step = 2,
    merged.tree.plot = TRUE,
    tree.node.label.split.character = NULL,
    tree.node.shape = 'circle',
    tree.node.size = 30,
    tree.node.text.size = 0.5,
    merged.tree.node.size.scale = 1.25,
    merged.tree.node.text.size.scale = 2.5,
    merged.tree.cell.frac.ci = FALSE,
    merged.tree.clone.as.branch = TRUE,
    mtcab.event.sep.char = ',',
    mtcab.branch.text.size = 1,
    mtcab.branch.width = 0.3,
    mtcab.node.size = 3,
    mtcab.node.label.size = 1,
    mtcab.node.text.size = 1.5,
    cell.plot = TRUE, num.cells = 100,
    cell.border.size = 0.25,
    cell.border.color = 'black',
    clone.grouping = 'horizontal',
    scale.monoclonal.cell.frac = TRUE,
    show.score = FALSE,
    cell.frac.ci = TRUE,
    disable.cell.frac = FALSE,
    out.dir = 'output.TB22446', out.format = 'pdf',
    overwrite.output = TRUE,
    width = 11, height = 7,
    panel.widths = c(3,4,2,4,4))

### fishy
f = generateFishplotInputs(results=y)
fishes = createFishPlotObjects(f)
#plot with fishplot
pdf(file = "fishplot.src172.mean.reagional.pdf", width=6, height=5)
for (i in 1:5){
fish = layoutClones(fishes[[i]])
fish = setCol(fish,f$clonevol.clone.colors)
fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
vlines=seq(1, length(vaf.col.names)), vlab=vaf.col.names, pad.left=0.5)
}
dev.off()



