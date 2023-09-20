#!/bin/sh
## run example GISTIC analysis ##

## output directory
echo --- creating output directory ---
basedir=`pwd`/example_results
mkdir -p $basedir 

echo --- running GISTIC ---
## input file definitions
#segfile=`pwd`/examplefiles/segmentationfile_9samples_NA.txt
segfile=/home/hmei/FreeC.txt
refgenefile=`pwd`/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
./gistic2 -b $basedir -seg $segfile -refgene $refgenefile
