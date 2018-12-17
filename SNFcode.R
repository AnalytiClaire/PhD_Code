##### SIMILARITY NETWORK FUSION #####

library(SNFtool)
library(gdata)

## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

#Read in network nodes
DEG_PPI <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes.txt")

#Extract PPI network genes from each dataset
C9 <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/C9uniquegene_samples.csv", row.names = 1)
C9 <- C9[,4:11]
C9 <- subset(C9, rownames(C9) %in% DEG_PPI)

VCP <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/VCPuniquegene_samples.csv", row.names = 1)
VCP <- VCP[,4:1]
VCP <- subset(VCP, rownames(VCP) %in% DEG_PPI)

FTLD <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/FTLDuniquegene_samples.csv", row.names = 1)
FTLD <- FTLD[,9:24]
FTLD <- subset(FTLD, rownames(FTLD) %in% DEG_PPI)

sALS <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/sALSuniquegene_samples.csv", row.names = 1)
sALS <- sALS[,4:10]
sALS <- subset(sALS, rownames(sALS) %in% DEG_PPI)

# PET <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/PET_results_keepfiltering.csv")
# rownames(PET) <- PET$hgnc_symbol
# PET <- PET[,19:35]
# PET <- subset(PET, rownames(PET) %in% DEG_PPI)
# 
# RAV <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/RAV_results_keepfiltering.csv")
# rownames(RAV) <- RAV$hgnc_symbol
# RAV <- RAV[,18:30]
# RAV <- subset(RAV, rownames(RAV) %in% DEG_PPI)

#Find the gene names that all datasets have in common
DEG_com <- Reduce(intersect, list(rownames(C9), rownames(VCP), rownames(FTLD),
                                  rownames(sALS)))

#Subset each dataset with these common names so they are all the same size
C9 <- subset(C9, rownames(C9) %in% DEG_com)
VCP <- subset(VCP, rownames(VCP) %in% DEG_com)
FTLD <- subset(FTLD, rownames(FTLD) %in% DEG_com)
sALS <- subset(sALS, rownames(sALS) %in% DEG_com)
# PET <- subset(PET, rownames(PET) %in% DEG_com)
# RAV <- subset(RAV, rownames(RAV) %in% DEG_com)

#Order by genename alphabetically
Data1 <- C9
Data1 <- Data1[order(rownames(Data1)),]
Data2 <- VCP
Data2 <- Data2[order(rownames(Data2)),]
Data3 <- FTLD
Data3 <- Data3[order(rownames(Data3)),]
Data4 <- sALS
Data4 <- Data4[order(rownames(Data4)),]
# Data5 <- PET
# Data5 <- Data5[order(rownames(Data5)),]
# Data6 <- RAV
# Data6 <- Data6[order(rownames(Data6)),]


## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows; if the data is discrete, we recommend the users to use "chiDist2"
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));
Dist4 = dist2(as.matrix(Data4),as.matrix(Data4));
# Dist5 = dist2(as.matrix(Data5),as.matrix(Data5));
# Dist6 = dist2(as.matrix(Data6),as.matrix(Data6));

#Returns an affinity matrix that represents the neighborhood graph of the data points.
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)
W3 = affinityMatrix(Dist3, K, alpha)
W4 = affinityMatrix(Dist4, K, alpha)
# W5 = affinityMatrix(Dist5, K, alpha)
# W6 = affinityMatrix(Dist6, K, alpha)

#Combine affinity matrices into fused network
W = SNF(list(W1,W2,W3,W4), K, T)

#name rows and columns
rownames(W) <- colnames(W) <- rownames(Data1)
# 
# ## spectral clustering
# C = 12 # number of clusters
# group = spectralClustering(W, C); 	# the final subtypes information
# colors <- c("red", "orange", "yellow", "green", "blue", "violet", 
#             "magenta", "brown", "grey","black", "teal", "white")
# col <- getColorsForGroups(group, colors)
# col <- data.frame(rownames(W), col)

#WGCNA clustering
library(WGCNA)
library(dynamicTreeCut)
dissTOM = TOMdist(as.matrix(W))
gsTree = hclust(as.dist(dissTOM), method = "average");

minClusterSize = 7
dynamicMods = cutreeDynamic(dendro = gsTree, distM = dissTOM, minClusterSize = minClusterSize,
                            cutHeight = "tree", deepSplit = 0)
sizeGrWindow(8,6)
plot(dynamicMods)
dynmods <- as.data.frame(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
colors <- table(dynamicColors)


geneInfo = data.frame(Gene = rownames(W),
                       module = dynamicColors)
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/SimilarityNetworkFusion/microarray/WGCNA_DS0_Min7")
# write.csv(W, "SNFoutput.csv")
write.csv(geneInfo, "geneInfo.csv", row.names = F)

# tdpgenes <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")
# tdpgenes <- tdpgenes$V1
# 
# braingenes <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv", header = T)
# allbrain <- braingenes$Gene.symbol
# 
# merge <- merge(geneInfo, braingenes, by.x = "Gene", by.y = "Gene.symbol")
# merge$braingene <- "TRUE"
# write.csv(merge, "braingene.csv", row.names = F, quote = F)

## Generate file for Cytoscape
#Remove doubled values
df <- W
df[lower.tri(df, diag = TRUE)] <- NA

#Turn into vector
df <- as.matrix(df)
dfvec <- unmatrix(df)
#Remove any NA values
dfvec <- na.omit(dfvec)
dfvec <- as.data.frame(dfvec)

#Split row names into two and assign to columns
library(stringr)
genenames <- as.data.frame(str_split_fixed(rownames(dfvec), ":", 2))
dfvec$gene1 <- genenames$V1
dfvec$gene2 <- genenames$V2

#reorder columns?me
dfvec <- dfvec[,c(2,3,1)]

#Take subset of edges if required
cyt <- subset(dfvec, dfvec > 0.003)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/SimilarityNetworkFusion/microarray/WGCNA_DS0_Min7")
for (i in 1:length(colors)){
  ref = names(colors[i])
  module = subset(geneInfo, module == ref)
  genelist <- module$Gene
  write.table(genelist, file = paste(names(colors[i]), "genelist.txt", sep = ""), row.names = F, quote = F, col.names = F)
  write.csv(module, file = paste(names(colors[i]), ".csv", sep = ""), row.names = F, quote = F)
}

ColList <- list()
for (i in 1:length(colors)){
  ref = colors[i]
  module = subset(geneInfo, module == ref)
  ColList[[i]] <- module$Gene
  names(ColList)[i] <- colors[i]
}

