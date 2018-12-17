#### Method for extracting genes from expression tables and taking the average across samples #####

library(matrixStats)

#Read in files
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")

A1 <- read.csv("C9rankeduniqueresult.csv")
A2 <- read.csv("CHrankeduniqueresult.csv")
A3 <- read.csv("sALSrankeduniqueresult.csv")
A4 <- read.csv("FTLDrankeduniqueresult.csv")
A5 <- read.csv("VCPrankeduniqueresult.csv")

#Read in genes of interest
setwd("/Users/clairegreen/Desktop/")
names <- read.table("TDP-43DEGs.txt")
names <- names$V1

#Take the subset of genes from the tables
subsetC9 <- subset(A5, A5$Gene.Symbol %in% names, drop = TRUE)
subsetC9 <-subsetC9[!duplicated(subsetC9[,15]),]
rownames(subsetC9) <- subsetC9$Gene.Symbol

#Take one condition (patients or controls)
subsetC9 <- subsetC9[,63:72]
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Expression_DEGonly")
write.csv(x = subsetC9, file = "sftld_DEG_expression.csv")

#Create empty column
subsetC9$vcp_median <- 0
#Convert to matrix so median will work
subsetC9 <- as.matrix(subsetC9)
#calculate median and add to new column
subsetC9[,ncol(subsetC9)] <- rowMedians(subsetC9) 

#Take only the rownames and the median value
result <- subsetC9[,ncol(subsetC9), drop = FALSE]

#subsetC9[,(ncol(subsetC9)-1)] <- rownames(subsetC9)
#subsetC9 <- subsetC9[,(ncol(subsetC9)-1):ncol(subsetC9)]

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/MedianGenes")
write.csv(result, file = "_mediansubset.csv")
