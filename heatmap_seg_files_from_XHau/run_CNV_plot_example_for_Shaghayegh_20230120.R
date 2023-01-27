

################################################################################

rm(list = ls())

##########

#hg19 constant:
#(modification needed if not hg19)

chrName <- c(1:22, "X", "Y")
nChr <- length(chrName)
chrLength <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566) #length of each chromosome
names(chrLength) <- chrName
#chrLength
#sum(chrLength)
#sum(chrLength[1:22])

##########

#Input data:
#Column: Subject_ID, Sample_ID, Chr, Start, End, Status.
#Each row represents one segment of one sample from one subject.
#The names of chromosome should be among c(1:22, "X", "Y")
#The types are "Amplification", "Deletion", "LOH", "Normal" here. Make any modificatin if needed.
#The orders of the rows don't matter.

dataSegment <- read.table("dataSegment_example_for_Shaghayegh_20230120.txt", header = TRUE, sep = "\t")
dataSegment$Chr <- as.character(dataSegment$Chr)
dataSegment <- dataSegment[!dataSegment$Chr %in% c("X", "Y"), ] #here we exclude the CNV on chromosome X and Y
#str(dataSegment)
#head(dataSegment)

nSampleBySubject <- rowSums(table(dataSegment$Subject_ID, dataSegment$Sample_ID) > 0) #number of samples from each subject
nSample <- sum(nSampleBySubject) #total number of samples

##########

#Generate data for plot:
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
colorSet <- c(Amplification = "red", Deletion = "blue", LOH = "yellow", Normal = "white")
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
par(cex.axis = 0.8, las = 1, mar = c(2, 2, 0, 0)) #parameters work well in current settings, modification needed in some cases

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
		tempXLeft <- 2.9e9 + c(0, cumsum(tempLength[c("Deletion", "LOH")])) / 2.9 * 0.4
		tempXRight <- 2.9e9 + cumsum(tempLength[c("Deletion", "LOH", "Amplification")]) / 2.9 * 0.4
		rect(tempXLeft, currentY - j + 0.15, tempXRight, currentY - j + 0.85, col = colorSet[c("Deletion", "LOH", "Amplification")], border = NA)
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
mtext(c("Chr", 1:22), 1, at = c(0, c(0, cumsum(chrLength[1:21])) + chrLength[1:22] / 2), cex = 0.9, line = -2)


#left subject IDs
tempNSampleBySubjectVec <- nSampleBySubject[names(dataPlotList)]
mtext(names(dataPlotList), 2, at = -c(0, cumsum(tempNSampleBySubjectVec)[-length(tempNSampleBySubjectVec)]) - tempNSampleBySubjectVec / 2, cex = 0.5, line = -5)

#topright axis
axis(3, seq(2.9e9, 3.3e9, 0.1e9), seq(0, 100, 25), pos = 1)
mtext("Proportion of genome (%)", 3, at = 3.1e9, line = -1, cex = 0.8)

#top legend
rect(c(0.5e9, 0.85e9, 1.5e9, 2e9), 1, c(0.55e9, 0.9e9, 1.55e9, 2.05e9), 3, pch = 22, col = colorSet[c("Deletion", "LOH", "Amplification", "Normal")])
text(c(0.55e9, 0.9e9, 1.55e9, 2.05e9), 2, c("Deletion", "Loss of heterozygosity", "Amplification", "Copy neutral"), pos = 4, cex = 1)

dev.off()

#####################################################################################



