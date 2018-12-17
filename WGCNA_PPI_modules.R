##### 
###WGCNA###
Genelist <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes.txt")

#load normalised Microarray datasets
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")
C9 <- read.csv("C9_unique.csv", row.names = 1)
C9gene <- subset(C9,rownames(C9) %in% Genelist )
sals <- read.csv("sals_unique.csv", row.names = 1)
salsgene <- subset(sals,rownames(sals) %in% Genelist)
FTLD <- read.csv("ftld_unique.csv", row.names = 1)
FTLDgene <- subset(FTLD,rownames(FTLD) %in% Genelist)
VCP <- read.csv("vcp_unique.csv", row.names = 1)
VCPgene <- subset(VCP,rownames(VCP) %in% Genelist)

#Load RNA-seq datasets
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
PET <- read.csv("PET_results_keepfiltering.csv", row.names = 36)
PETgene <- subset(PET,rownames(PET) %in% Genelist)
RAV <- read.csv("RAV_results_keepfiltering.csv", row.names = 31)
RAVgene <- subset(RAV, rownames(RAV) %in% Genelist)

C9pat <- C9gene[,11:18]
salspat <- salsgene[,11:17]
sFTLDpat <- FTLDgene[,22:31]
GRN_FTLDpat <- FTLDgene[16:21]
VCPpat <- VCPgene[,11:17]
C9_PETpat <- PETgene[,1:7]
RAVpat <- RAVgene[,18:30]

library(WGCNA)
options(stringsAsFactors = FALSE);

Analysis.name <- "VCP"


Exp1 <- VCPpat

Exp <- t(Exp1)

###PATIENT ANALYSIS###

###Choosing soft threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(Exp, powerVector = powers, verbose = 5)
# Plot the results:
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.50,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##SOFT THRESHOLD VALUE SELECTED##

# Exp <- data.matrix(Exp) #csv files contain character matrices, the following code requires numeric

##One-step network construction and module detection
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_WGCNA/VCP/")
net = blockwiseModules(Exp, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)



dyntree <- dynamicTreeCut::cutreeDynamicTree(net$dendrograms[[1]], maxTreeHeight = 1, deepSplit = TRUE, minModuleSize = 50)




table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "-networkConstruction-auto.RData")


#find out what the IDs are of the genes that are contained within a module. 

colnames(Exp)[moduleColors=='green']


geneInfo0 = data.frame(ProbeID = rownames(Exp1),
                       moduleColor = moduleColors)
rownames(geneInfo0) <- NULL



###Exporting to Cytoscape

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
x <- geneInfo0$moduleColor
x <- unique(x)

# intModules <- c("green", "purple", "brown", "greenyellow", "royalblue")
TOM = TOMsimilarityFromExpr(Exp, power = 10);

for (i in 1:length(x)){
  intModules <- x[[i]]

  probes = colnames(Exp)
  inModule = is.finite(match(moduleColors, intModules));
  modProbes = probes[inModule];
  modGenes = geneInfo0$geneSymbol [match(modProbes, geneInfo0$ProbeID)];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(intModules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(intModules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]);
}



# C9out <- split(geneInfo0, f = geneInfo0$moduleColor)
C9list <- C9out
C9_colours <- lapply(C9list, '[[', 1)

# sALSout <- split(geneInfo0 , f = geneInfo0$moduleColor)
sALSlist <- sALSout
sALS_colours <- lapply(sALSlist, '[[', 1)

# sFTLDout<- split(geneInfo0 , f = geneInfo0$moduleColor)
sFTLDlist <- sFTLDout
sFTLD_colours <- lapply(sFTLDlist, '[[', 1)

# GRN_FTLDout <-split(geneInfo0 , f = geneInfo0$moduleColor)
GRN_FTLDlist <- GRN_FTLDout
GRN_FTLD_colours <- lapply(GRN_FTLDlist, '[[', 1)

# VCPout <-split(geneInfo0 , f = geneInfo0$moduleColor)
VCPlist <- VCPout
VCP_colours <- lapply(VCPlist, '[[', 1)

list <- C9out

for (i in 1:length(list)){
  x <- as.data.frame(list[i])
}

df <- data.frame(matrix(unlist(C9out), byrow=T))








#####
### Find overlap of modules using Jaccard coefficient
sub_list1 <- sFTLD_colours #make sure this is the bigger list
sub_list2 <- sALS_colours #make sure this is the smaller list

#Make function for Jaccard coefficient
jaccard<-function(A,B){
  jc=length(intersect(A,B))/length(union(A,B))
  return(jc)
}

#Create empty matrix 
jc_mat <- matrix(0, length(sub_list1), length(sub_list2))

#Compare modules in sub_list1 and sub_list2
for (i in 1:length(sub_list2)) {
  list <- sub_list2[[i]]
  for (j in 1:length(sub_list1)){
    jc_mat[j,i] <- jaccard(sub_list1[[j]], list)
  }
}

rownames(jc_mat) <- names(sub_list1)
colnames(jc_mat) <- names(sub_list2)

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_WGCNA/")
write.csv(jc_mat, "row_sFTLD_col_sALS.csv")


