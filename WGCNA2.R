library(WGCNA)
options(stringsAsFactors = FALSE);
### orf72 ###
# Display the current working directory
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/WGCNA/C9orf72/")

Analysis.name <- "C9"

#Read in desired genes
Results <- read.csv ("C9result.csv", header=TRUE, row.names = 1) 
# #Results <- Results[1:8000,]
# #gene expression analysis to find criteria
# 
# ID <- as.data.frame(Results$Probe.Set.ID)
# colnames(ID)[1] <- "ProbeID"


Exp <- Results
  #merge(ID, RawExp, by.x="ProbeID", by.y="Probe.Set.ID") #merge raw expression data with accepted genes
# rownames(Exp) <- Exp[,15] #make probeset IDs row names
#colnames(Exp) <- colnames(RawExp) #make file names column names
# Exp[,1] <- NULL  #remove ID column
# Exp <- Exp[,51:58] #Take only patients

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
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/WGCNA/C9orf72")
net = blockwiseModules(Exp, power = 4,
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
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "-networkConstruction-auto.RData")


### RELATING MODULES TO EXTERNAL INFORMATION AND IDENTIFYING IMPORTANT GENES ###
#lnames = load(file = "-networkConstruction-auto.RData")

# traitdata <- read.csv(file = "Pheno.csv")
# Samples = rownames(Exp)
# traitRows = match(Samples, traitdata$SampleName)
# datTraits = traitdata[traitRows, -1]
# rownames(datTraits) = datTraits[traitRows, 1]
# Traitpheno <- datTraits[,2]
# 
# # Define numbers of genes and samples
# nGenes = ncol(Exp)
# nSamples = nrow(Exp)
# # Recalculate MEs with color labels
# MEs0 = moduleEigengenes(Exp, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)
# moduleTraitCor = cor(MEs, Traitpheno, use = "p")
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 
# 
# sizeGrWindow(10,6)
# # Will display correlations and their p-values
# textMatrix = paste(signif(moduleTraitCor, 2), "(", signif(moduleTraitPvalue, 1), ")", sep = " ");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3))
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = "Disease Phenotype",
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = greenWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.9,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# 
# 
# ### Gene relationship to trait and important modules: Gene Significance and Module Membership ##
# 
# # Define variable weight containing the weight column of datTrait
# pheno = as.data.frame(datTraits$PhenoLabel)
# names(pheno) = "Disease Phenotype"
# # names (colors) of the modules
# modNames = substring(names(MEs), 3)
# geneModuleMembership = as.data.frame(cor(Exp, MEs, use = "p"))
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# names(geneModuleMembership) = paste("MM", modNames, sep="")
# names(MMPvalue) = paste("p.MM", modNames, sep="")
# geneTraitSignificance = as.data.frame(cor(Exp, pheno, use = "p"))
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# names(geneTraitSignificance) = paste("GS.", names(pheno), sep="")
# names(GSPvalue) = paste("p.GS.", names(pheno), sep="")
# 
# module = "brown"
# column = match(module, modNames)
# moduleGenes = moduleColors==module
# sizeGrWindow(7, 7)
# par(mfrow = c(1,1))
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for disease phenotype",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#find out what the IDs are of the genes that are contained within a module. 

colnames(Exp)[moduleColors=='green']


geneInfo0 = data.frame(ProbeID = Results$Probe.Set.ID,
                       geneSymbol = Results$Gene.Symbol,
                       moduleColor = moduleColors)
rownames(geneInfo0) <- NULL

modOrder = order(-abs(cor(MEs, pheno, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")


### Interfacing network analysis with other data such as functional annotation and gene ontology ###

## Output gene lists for use with online software and services

Entrez <- cbind(Results$Probe.Set.ID, Results$Entrez.Gene) 
colnames(Entrez)[1] <- "ProbeID"
colnames(Entrez)[2] <- "EntrezID"

geneInfo <- merge(geneInfo, Entrez, by.x="ProbeID", by.y="ProbeID")

AllIDs <- geneInfo$geneSymbol
# # $ Choose interesting modules
# intModules = c("royalblue", "brown", "green", "purple", "greenyellow")
# for (module in intModules)
# {
#   # Select module probes
#   modGenes = (moduleColors==module)
#   # Get their entrez ID codes
#   modLLIDs = IDs[modGenes];
#   # Write them into a file
#   fileName = paste("IDs-", module, ".txt", sep="");
#   write.table(as.data.frame(modLLIDs), file = fileName,
#               row.names = FALSE, col.names = FALSE)
# }
# # As background in the enrichment analysis, we will use all probes in the analysis.
# fileName = paste("IDs-all.txt", sep="");
# write.table(as.data.frame(IDs), file = fileName,
#             row.names = FALSE, col.names = FALSE)

#Enrichment analysis directly within R
source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
library(org.Hs.eg.db)

GOenr = GOenrichmentAnalysis(moduleColors, AllIDs, ontologies = c("BP"), organism = "human", nBestP = 10)

tab = GOenr$bestPTerms[[4]]$enrichment

write.table(tab, file = "GOEnrichmentTableBP.csv", sep = ",", quote = TRUE, row.names = FALSE)


###Exporting to Cytoscape

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
intModules <- c("green", "purple", "brown", "greenyellow", "royalblue")

intModules <- "brown"

TOM = TOMsimilarityFromExpr(Exp, power = 6);
probes = colnames(Exp)
inModule = is.finite(match(moduleColors, intModules));
modProbes = probes[inModule];
modGenes = geneInfo$geneSymbol [match(modProbes, geneInfo$ProbeID)];
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



# dissTOM = 1-TOMsimilarityFromExpr(Exp, power = 6);
# # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
# plotTOM = dissTOM^7;
# # Set diagonal to NA for a nicer plot
# diag(plotTOM) = NA;
# 
# # Call the plot function
# sizeGrWindow(9,9)
# TOMplot(plotTOM, moduleColors, main = "Network heatmap plot, all genes")
# 
# # Recalculate module eigengenes
# MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# # Isolate weight from the clinical traits
# weight = as.data.frame(datTraits$weight_g);
# names(weight) = "weight"
# # Add the weight to existing module eigengenes
# MET = orderMEs(cbind(MEs, weight))
# # Plot the relationships among the eigengenes and the trait
# sizeGrWindow(5,7.5);
# par(cex = 0.9)
# plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
#                       = 90)