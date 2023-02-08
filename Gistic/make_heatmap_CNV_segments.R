
# Title: marge CNV segments for each subject ID 
# Original Author:  Xing Hau
# Contributors:    Shaghayegh Soudi 
# Date: January 2023

### load required libraries
#library(bedtoolsr)
#library(GenomicRanges)
library(data.table)
library(VariantAnnotation)
meta<-read.table(file = "/scratch/users/shsoudi/metadata/WES_updated_metadat_02022023.txt", header = TRUE)
meta$RT_status[meta$sequenceofsamplevRT=="beforeRT"]<-"noRT"
meta$RT_status[meta$sequenceofsamplevRT=="nopreopRT"]<-"noRT"
meta$RT_status[meta$sequenceofsamplevRT=="afterRT"]<-"afterRT"

seg_files<-list.files("00-inputs-segs", pattern = "*.titan.ichor.segfull.txt", full.names = TRUE)

attackTitan1 <- function(x) {

      data<-read.delim(x,  header=TRUE, sep = "\t")
      data$sample_id<-gsub("-.*$","",data$Sample)
      data<-merge(data,meta, by.x = "sample_id", by.y = "sampleid")
      data$uniq_id<-gsub("_.*$","",data$sample_id)
       data$uniq_id2<-paste(data$uniq_id,data$RT_status, sep = "_")
      #data$major_minor_copy<-paste(data$MajorCN,data$MinorCN, sep = "-")
      return(data)
    }

titan_segs <- lapply(seg_files, attackTitan1)
titan_segs <- do.call("rbind", titan_segs)
patients_rt<-unique(titan_segs$uniq_id2)

dataSegment<-titan_segs[,c("uniq_id2","sample_id","Chromosome", "Start", "End", "Corrected_Call")]
colnames(dataSegment)<-c("Subject_ID","Sample_ID","Chr","Start","End","Corrected_Call")

#### run main part to merge segments for each subect ID
titan_segs<-dataSegment[dataSegment$Corrected_Call !="Fill_In",]
titan_segs$type[titan_segs$Corrected_Call=="NEUT"]<-"Normal"
titan_segs$type[titan_segs$Corrected_Call== "AMP"| titan_segs$Corrected_Call== "GAIN" | titan_segs$Corrected_Call== "HLAMP"]<-"Amplification"
titan_segs$type[titan_segs$Corrected_Call== "HETD"| titan_segs$Corrected_Call== "HOMD"]<-"Deletion"

names(titan_segs)[7]<-"Status"
titan_segs<-titan_segs[,c("Subject_ID","Sample_ID","Chr","Start" ,"End","Status")]

titan_segs$Chr<-gsub("chr","",titan_segs$Chr)
titan_segs$Chr<-gsub("X","23",titan_segs$Chr)

write.table(titan_segs, file = "dataSegment_for_segplot_02062023.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


dataSegment<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_02062023/dataSegment_for_segplot_02062023.txt", header = TRUE)
chrName <- c(1:23)
nChr <- length(chrName)
chrLength <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560) #length of each chromosome
names(chrLength) <- chrName
#chrLength
#sum(chrLength)
#sum(chrLength[1:22])

#dataSegment <- dataSegment[!dataSegment$Chr %in% c("X", "Y"), ] #here we exclude the CNV on chromosome X and Y
nSampleBySubject <- rowSums(table(dataSegment$Subject_ID, dataSegment$Sample_ID) > 0) #number of samples from each subject
nSample <- sum(nSampleBySubject) #total number of samples
dataPlot <- dataSegment

#Locate the position of the each segment on the whole genome:
tempLengthVec <- cumsum(c(0, chrLength[-nChr])); names(tempLengthVec) <- chrName
dataPlot$xLeft <- dataPlot$Start + tempLengthVec[dataPlot$Chr]
dataPlot$xRight <- dataPlot$End + tempLengthVec[dataPlot$Chr]
dataPlot$Length <- dataPlot$xRight - dataPlot$xLeft + 1

##########

#Set the color for each type of the segments (make sure the names match those in column 'Status' of dataSegment):
#(change the colors if neeeded)
colorSet <- c(Amplification = "red", Deletion = "blue",  Normal = "white")
dataPlot$color <- colorSet[dataPlot$Status]

##########

dataPlot <- dataPlot[order(dataPlot$Start), ]
dataPlot <- dataPlot[order(dataPlot$Chr), ]
dataPlot <- dataPlot[order(dataPlot$Sample_ID), ]
dataPlot <- dataPlot[order(dataPlot$Subject_ID), ]
#str(dataPlot)
#head(dataPlot)

##########

#Each element of the list 'dataPlotList' is the segments of the samples from one subject.
#The plot will follow the order of the list (from top to bottom).
#Here we use the order by number of samples from each subject (change the order of the subjcts if necessary).

dataPlotList <- split(dataPlot, dataPlot$Subject_ID)
dataPlotList <- dataPlotList[names(sort(nSampleBySubject, decreasing = TRUE))]
#length(dataPlotList)
#names(dataPlotList)

#####################################################################################

#Generate the figure based on dataPlotList:
#The main steps are: (1) generate the blank background; (2) adding rectangles to the correspoding locations; (3) adding labels/legends.

png("CNV_plot_example_data_for_Shaghayegh_02062023.png", width = 3600, height = 4800, res = 300) #change the file name, size or resulotion of the figure if necessary
par(cex.axis = 0.5, las = 1, mar = c(2, 2, 2, 0)) #parameters work well in current settings, modification needed in some cases

#the total length of the chorosome (exclude X and Y) is around 2.9e9
#the right bars locate from 2.9e9 to 3.3e9 on x-axis
#the left lables of sample/subject IDs located from -0.2e9 to 0 on x-axis
#the i-th sample locates from -(i - 1) to -i on y-axis
plot(NA, xlim = c(-0.2e9,3.3e9), ylim = c(-nSample, 3), axes = FALSE, xlab = "", ylab = "")

currentY <- 0

for(i in 1:length(dataPlotList)){
	tempSubjectID <- names(dataPlotList)[i]
	tempDataSegments <- dataPlotList[[i]]
	tempSampleID <- sort(unique(tempDataSegments$Sample_ID))
	tempNSample <- length(tempSampleID)
	tempSampleIDFactor <- factor(tempDataSegments$Sample_ID, levels = tempSampleID)

#	adding CNV segments:
	tempDataSegments$yTop <- currentY - as.integer(tempSampleIDFactor) + 1
	tempDataSegments$yBottom <- currentY - as.integer(tempSampleIDFactor)
	rect(tempDataSegments$xLeft, tempDataSegments$yBottom, tempDataSegments$xRight, tempDataSegments$yTop, col = tempDataSegments$color, border = NA)

#	adding right bars:
	tempAggregate <- aggregate(tempDataSegments$Length, list(tempDataSegments$Sample_ID, tempDataSegments$Status), sum, drop = FALSE)
	for(j in 1:tempNSample){
		tempLength <- tempAggregate$x[tempAggregate$Group.1 == tempSampleID[j]]
		names(tempLength) <- tempAggregate$Group.2[tempAggregate$Group.1 == tempSampleID[j]]
		tempXLeft <- 2.9e9 + c(0, cumsum(tempLength[c("Deletion")])) / 2.9 * 0.4
		tempXRight <- 2.9e9 + cumsum(tempLength[c("Deletion", "Amplification")]) / 2.9 * 0.4
		rect(tempXLeft, currentY - j + 0.15, tempXRight, currentY - j + 0.85, col = colorSet[c("Deletion", "Amplification")], border = NA)
	}

	currentY <- currentY - tempNSample
}

#if(currentY != -nSample)stop("Something wrong!")

tempX <- sum(chrLength[1:22])
for(i in -(0:nSample))lines(c(0, tempX), c(i, i), col = "gray", lwd = 0.1) #thin horizontal line
tempYVec <- c(0, -cumsum(nSampleBySubject[names(dataPlotList)]))
for(i in 1:length(tempYVec))lines(c(-0.2e9,3.3e9), c(tempYVec[i], tempYVec[i]), col = "black", lwd = 1) #thick horizontal line
lines(c(0, 0), c(0, -nSample), col = "black", lwd = 1) #left vertical line
lines(c(tempX, tempX), c(0, -nSample), col = "black", lwd = 1) #right vertical line
for(i in 1:21)lines(c(cumsum(chrLength)[i], cumsum(chrLength)[i]), c(0, -nSample), col = "black", lwd = 0.5) #vertical lines seperate chromosome

#bottom chromosome name
mtext(c("Chr", 1:22), 1, at = c(2, c(0, cumsum(chrLength[1:21])) + chrLength[1:22] / 2), cex = 0.4, line = -2)


#left subject IDs
tempNSampleBySubjectVec <- nSampleBySubject[names(dataPlotList)]
mtext(names(dataPlotList), 2, at = -c(2, cumsum(tempNSampleBySubjectVec)[-length(tempNSampleBySubjectVec)]) - tempNSampleBySubjectVec / 2, cex = 0.5, line = -5)

#topright axis
axis(3, seq(2.9e9, 3.3e9, 0.1e9), seq(0, 100, 25), pos = 1)
mtext("Proportion of genome (%)", 3, at = 3.1e9, line = -1, cex = 0.4)

#top legend
rect(c(0.5e9, 0.85e9, 1.5e9, 2e9), 1, c(0.55e9, 0.9e9, 1.55e9, 2.05e9), 3, pch = 22, col = colorSet[c("Deletion", "LOH", "Amplification", "Normal")])
text(c(0.55e9, 0.9e9, 1.55e9, 2.05e9), 2, c("Deletion", "Amplification", "Copy neutral"), pos = 4, cex = 0.5)

dev.off()

#####################################################################################




########
#########
##########


meta<-read.table(file = "/scratch/users/shsoudi/metadata/metadata.txt", header = FALSE)
seg_files<-list.files("00-inputs-segs", pattern = "*.segs.txt", full.names = TRUE)

attackTitan1 <- function(x) {

      data<-read.delim(x,  header=TRUE, sep = "\t")
      data$sample_id<-gsub("-.*$","",data$Sample)
      data<-merge(data,meta, by.x = "sample_id", by.y = "V1")
      data$uniq_id<-gsub("_.*$","",data$sample_id)
      data$uniq_id2<-paste(data$uniq_id,data$V5, sep = "_")
      #data$major_minor_copy<-paste(data$MajorCN,data$MinorCN, sep = "-")
      return(data)
    }

titan_segs <- lapply(seg_files, attackTitan1)
titan_segs <- do.call("rbind", titan_segs)
patients_rt<-unique(titan_segs$uniq_id2)


dataSegment<-titan_segs[,c("uniq_id2","sample_id","Chromosome", "Start_Position.bp.", "End_Position.bp.", "TITAN_call")]
colnames(dataSegment)<-c("Subject_ID","Sample_ID","Chr","Start","End","Status")


dataSegment$Status<-gsub("NLOH","Normal",dataSegment$Status)
dataSegment$Status<-gsub("HET","Normal",dataSegment$Status)
dataSegment$Status<-gsub("DLOH","LOH",dataSegment$Status)
dataSegment$Status<-gsub("ALOH","LOH",dataSegment$Status)
dataSegment$Status<-gsub("GAIN","Amplification",dataSegment$Status)
dataSegment$Status<-gsub("UBCNA","Amplification",dataSegment$Status)
dataSegment$Status<-gsub("HOMD","Deletion",dataSegment$Status)
dataSegment$Status<-gsub("BCNA","Amplification",dataSegment$Status)
dataSegment$Status<-gsub("ASCNA","Amplification",dataSegment$Status)

dataSegment$Chr<-gsub("chr","",dataSegment$Chr)
write.table(dataSegment, file = "dataSegment.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


dataSegment<-read.table(file = "~/Desktop/dataSegment.txt", header = TRUE)


dataSegment$Status<-gsub("LOH","Deletion",dataSegment$Status)


chrName <- c(1:22, "X")
nChr <- length(chrName)
chrLength <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566) #length of each chromosome
names(chrLength) <- chrName
#chrLength
#sum(chrLength)
#sum(chrLength[1:22])

dataSegment <- dataSegment[!dataSegment$Chr %in% c("X", "Y"), ] #here we exclude the CNV on chromosome X and Y
nSampleBySubject <- rowSums(table(dataSegment$Subject_ID, dataSegment$Sample_ID) > 0) #number of samples from each subject
nSample <- sum(nSampleBySubject) #total number of samples


dataPlot <- dataSegment


##########

#Locate the position of the each segment on the whole genome:
tempLengthVec <- cumsum(c(0, chrLength[-nChr])); names(tempLengthVec) <- chrName
dataPlot$xLeft <- dataPlot$Start + tempLengthVec[dataPlot$Chr]
dataPlot$xRight <- dataPlot$End + tempLengthVec[dataPlot$Chr]
dataPlot$Length <- dataPlot$xRight - dataPlot$xLeft + 1

##########

#Set the color for each type of the segments (make sure the names match those in column 'Status' of dataSegment):
#(change the colors if neeeded)
colorSet <- c(Amplification = "red", Deletion = "blue",  Normal = "white")
dataPlot$color <- colorSet[dataPlot$Status]

##########

dataPlot <- dataPlot[order(dataPlot$Start), ]
dataPlot <- dataPlot[order(dataPlot$Chr), ]
dataPlot <- dataPlot[order(dataPlot$Sample_ID), ]
dataPlot <- dataPlot[order(dataPlot$Subject_ID), ]
#str(dataPlot)
#head(dataPlot)

##########

#Each element of the list 'dataPlotList' is the segments of the samples from one subject.
#The plot will follow the order of the list (from top to bottom).
#Here we use the order by number of samples from each subject (change the order of the subjcts if necessary).

dataPlotList <- split(dataPlot, dataPlot$Subject_ID)
dataPlotList <- dataPlotList[names(sort(nSampleBySubject, decreasing = TRUE))]
#length(dataPlotList)
#names(dataPlotList)

#####################################################################################

#Generate the figure based on dataPlotList:
#The main steps are: (1) generate the blank background; (2) adding rectangles to the correspoding locations; (3) adding labels/legends.

png("CNV_plot_example_data_for_Shaghayegh_20230120.png", width = 3600, height = 4800, res = 300) #change the file name, size or resulotion of the figure if necessary
par(cex.axis = 0.5, las = 1, mar = c(2, 2, 2, 0)) #parameters work well in current settings, modification needed in some cases

#the total length of the chorosome (exclude X and Y) is around 2.9e9
#the right bars locate from 2.9e9 to 3.3e9 on x-axis
#the left lables of sample/subject IDs located from -0.2e9 to 0 on x-axis
#the i-th sample locates from -(i - 1) to -i on y-axis
plot(NA, xlim = c(-0.2e9,3.3e9), ylim = c(-nSample, 3), axes = FALSE, xlab = "", ylab = "")

currentY <- 0

for(i in 1:length(dataPlotList)){
	tempSubjectID <- names(dataPlotList)[i]
	tempDataSegments <- dataPlotList[[i]]
	tempSampleID <- sort(unique(tempDataSegments$Sample_ID))
	tempNSample <- length(tempSampleID)
	tempSampleIDFactor <- factor(tempDataSegments$Sample_ID, levels = tempSampleID)

#	adding CNV segments:
	tempDataSegments$yTop <- currentY - as.integer(tempSampleIDFactor) + 1
	tempDataSegments$yBottom <- currentY - as.integer(tempSampleIDFactor)
	rect(tempDataSegments$xLeft, tempDataSegments$yBottom, tempDataSegments$xRight, tempDataSegments$yTop, col = tempDataSegments$color, border = NA)

#	adding right bars:
	tempAggregate <- aggregate(tempDataSegments$Length, list(tempDataSegments$Sample_ID, tempDataSegments$Status), sum, drop = FALSE)
	for(j in 1:tempNSample){
		tempLength <- tempAggregate$x[tempAggregate$Group.1 == tempSampleID[j]]
		names(tempLength) <- tempAggregate$Group.2[tempAggregate$Group.1 == tempSampleID[j]]
		tempXLeft <- 2.9e9 + c(0, cumsum(tempLength[c("Deletion")])) / 2.9 * 0.4
		tempXRight <- 2.9e9 + cumsum(tempLength[c("Deletion", "Amplification")]) / 2.9 * 0.4
		rect(tempXLeft, currentY - j + 0.15, tempXRight, currentY - j + 0.85, col = colorSet[c("Deletion", "Amplification")], border = NA)
	}

	currentY <- currentY - tempNSample
}

#if(currentY != -nSample)stop("Something wrong!")

tempX <- sum(chrLength[1:22])
for(i in -(0:nSample))lines(c(0, tempX), c(i, i), col = "gray", lwd = 0.1) #thin horizontal line
tempYVec <- c(0, -cumsum(nSampleBySubject[names(dataPlotList)]))
for(i in 1:length(tempYVec))lines(c(-0.2e9,3.3e9), c(tempYVec[i], tempYVec[i]), col = "black", lwd = 1) #thick horizontal line
lines(c(0, 0), c(0, -nSample), col = "black", lwd = 1) #left vertical line
lines(c(tempX, tempX), c(0, -nSample), col = "black", lwd = 1) #right vertical line
for(i in 1:21)lines(c(cumsum(chrLength)[i], cumsum(chrLength)[i]), c(0, -nSample), col = "black", lwd = 0.5) #vertical lines seperate chromosome

#bottom chromosome name
mtext(c("Chr", 1:22), 1, at = c(2, c(0, cumsum(chrLength[1:21])) + chrLength[1:22] / 2), cex = 0.4, line = -2)


#left subject IDs
tempNSampleBySubjectVec <- nSampleBySubject[names(dataPlotList)]
mtext(names(dataPlotList), 2, at = -c(2, cumsum(tempNSampleBySubjectVec)[-length(tempNSampleBySubjectVec)]) - tempNSampleBySubjectVec / 2, cex = 0.5, line = -5)

#topright axis
axis(3, seq(2.9e9, 3.3e9, 0.1e9), seq(0, 100, 25), pos = 1)
mtext("Proportion of genome (%)", 3, at = 3.1e9, line = -1, cex = 0.4)

#top legend
rect(c(0.5e9, 0.85e9, 1.5e9, 2e9), 1, c(0.55e9, 0.9e9, 1.55e9, 2.05e9), 3, pch = 22, col = colorSet[c("Deletion", "LOH", "Amplification", "Normal")])
text(c(0.55e9, 0.9e9, 1.55e9, 2.05e9), 2, c("Deletion", "Loss of heterozygosity", "Amplification", "Copy neutral"), pos = 4, cex = 0.5)

#dev.off()

#####################################################################################





