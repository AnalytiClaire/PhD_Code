library(WGCNA)
options(stringsAsFactors = FALSE);
### orf72 ###
# Display the current working directory
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/")

Analysis.name <- "VCP"

#Read in desired genes
Results <- read.csv ("VCPrankeduniqueresult.csv", header=TRUE) #Taking only the genes we deemed acceptable through criteria. See details on
Results <- Results[1:8000,]

Exp <- Results
rownames(Exp) <- Exp[,1] #make probeset IDs row names
#colnames(Exp) <- colnames(RawExp) #make file names column names
Exp[,1] <- NULL  #remove ID column
Exp <- Exp[,51:57] #Take only patients

# Pat <- Exp[,1:8]
# Con <- Exp[,9:11]
# 
# Pat <- t(Exp)
# Con <- t(Exp)
Exp <- t(Exp)

###PATIENT ANALYSIS###

###Choosing soft threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(Exp, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
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
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/")
net = blockwiseModules(Exp, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
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
# save(MEs, moduleLabels, moduleColors, geneTree,
#      file = "-networkConstruction-auto.RData")

geneInfo0 = data.frame(ProbeID = Results$Probe.Set.ID,
                       geneSymbol = Results$Gene.Symbol,
                       moduleColor = moduleColors)
rownames(geneInfo0) <- NULL

geneOrder = order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "VCPgeneInfo.csv")
