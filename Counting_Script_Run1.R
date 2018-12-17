#!/usr/bin/env Rscript

#Requires Rsamtools, GenomicFeatures and GenomicAlignments
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")

setwd("/home/elaine/Fastq/Run1/")
#Load sample table to obtain sample names and bam file information
csvfile <- "sample_table_TDP_Run1.csv"
sampleTable <- read.csv(csvfile,row.names=1)
filenames <- paste0(sampleTable$Run,sampleTable$SampleName, ".filtered_u.s.bam")

#Load GRanges datasets for which we get counts
load("/home/jenny/TDP-43_SortRNAseq/150727_Genes_ExonIntron_SharedRegionsRemoved.rda")

#Create count matrices to store counts into
counts.g <- matrix(nrow=length(genes1),ncol=length(filenames))
colnames(counts.g) <- sampleTable$SampleName
rownames(counts.g) <- names(genes1)

read.stats <- matrix(nrow=7,ncol=length(filenames))
rownames(read.stats) <- c("Total","MappingToGenes","MappingWithinGenes","Non-overlapping_Genic","MappingToEIBG","MappingWithinEIBG","GenicReadsMappingtoEIBG")
colnames(read.stats) <- sampleTable$SampleName

for(i in 1:length(filenames)) {	
	print(filenames[i])
	bam <- readGAlignmentPairs(filenames[i], use.names=TRUE)
	read.stats[1,i] <- length(bam)
	ov <- findOverlaps(bam,genes1, ignore.strand=TRUE) #finds reads that overlap with genes and files them such that if a read maps to 2 genes, it will be present twice (1 for each gene)
	read.stats[2,i] <- length(unique(queryHits(ov))) #includes reads that overlap with multiple genes
	ov.g <- findOverlaps(bam, genes1, type="within",ignore.strand=TRUE)
	read.stats[3,i] <- length(unique(queryHits(ov.g))) #excludes reads that overlap with genic and intergenic space
	reads_to_keep <- which(countQueryHits(ov.g)==1L)
	read.stats[4,i] <- length(reads_to_keep) #includes reads that only maps to one gene
	ov.g <- ov.g[queryHits(ov.g) %in% reads_to_keep]
	counts.g[,i] <- countSubjectHits(ov.g)
	
}

#Save count matrices and statistics
save(counts.g,read.stats,file="150809_Counts_SharedRegionsRemoved_Statistics_Run1.rda")