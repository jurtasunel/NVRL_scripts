#!/bin/bash
### This script aligns fastq files with bowtie2 and produces a text file with unmapped reads, a text file with depth information, ggplot of the depth and a csv file with low depth positions.
### This script requires an indexed reference genome and the Rscript plot_depth.R.

# Create a variable to store the reference path. Requires a previously created reference with bowtie2-build.
reference="/path/to/reference"
# Get the name of the fastq file as an argument from terminal.
rawfastq=$1

#COMANDS
#Indexing
#bowtie2-build '$param.ref' sarscov
# Align the input file to the reference with bowtie2 and write out a txt documment with unnaligned reads.
echo Aligning $1 with bowtie2...
bowtie2 --un ${rawfastq}_not_aligned.txt -x ${reference} -q ${rawfastq} -S ${rawfastq}.sam
# Convert the sam to bam file and remove sam file.
echo Converting sam to bam file...
samtools view -b ${rawfastq}.sam > ${rawfastq}.bam
rm ${rawfastq}.sam
# Sort the sam file and remove unsorted bam.
echo Sorting bam file...
samtools sort ${rawfastq}.bam -o ${rawfastq}_srtd.bam
rm ${rawfastq}.bam
# Index the sorted bam file.
echo Indexing bam file...
samtools index ${rawfastq}_srtd.bam
# Get depth information to text file and remove bam files.
echo Getting depth from bam file...
samtools depth -a -H ${rawfastq}_srtd.bam -o ${rawfastq}_depth.txt
rm ${rawfastq}_srtd.bam
rm ${rawfastq}_srtd.bam.bai
# Call Rscript for plotting.
Rscript plot_depth.R `pwd`/${rawfastq}_depth.txt



