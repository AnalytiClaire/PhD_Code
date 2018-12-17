setwd(dir="/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")

library(WGCNA)
library(impute)
library(dynamicTreeCut)
library(qvalue)
library(flashClust)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(Hmisc)

C9Expr <- read.csv(file = "C9rankeduniqueresult.csv")
rownames(C9Expr) <- C9Expr$Probe.Set.ID
datExprA1 <- C9Expr[,52:59] 

CHExpr <- read.csv(file = "CHrankeduniqueresult.csv")
rownames(CHExpr) <- CHExpr$Probe.Set.ID
datExprA2 <- CHExpr[,55:57] 

sALSExpr <- read.csv(file = "sALSrankeduniqueresult.csv")
rownames(sALSExpr) <- sALSExpr$Probe.Set.ID
datExprA3 <- sALSExpr[,52:58]

VCPExpr <- read.csv(file = "VCPrankeduniqueresult.csv")
rownames(VCPExpr) <- VCPExpr$Probe.Set.ID
datExprA4 <- VCPExpr[,52:58]

PETExpr <- read.csv(file = "PETEXPRSrankeduniqueresult.csv")
rownames(PETExpr) <- PETExpr$hgnc_symbol
datExprA5g <- PETExpr[,9:26]

RAVExpr <- read.csv(file = "RAVEXPRSrankeduniqueresult.csv")
rownames(RAVExpr) <- RAVExpr$hgnc_symbol
datExprA6g <- RAVExpr[,17:28]

FTLDExpr <- read.csv(file = "FTLDrankeduniqueresult.csv")
rownames(FTLDExpr) <- FTLDExpr$Probe.Set.ID
datExprA7 <- FTLDExpr[,57:72]



# CHExpr <- read.csv(file = "CHrankeduniqueresult.csv")
# rownames(CHExpr) <- CHExpr$Probe.Set.ID
# datExprB1 <- CHExpr[,55:57] #C9 patients
# datExprB2 <- CHExpr[,49:54] #C9 controls

IDs <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/C9orf72/C9ID geneIDs.csv")
IDs <- C9Expr[,c("Probe.Set.ID", "Gene.Symbol")]
genes <- IDs$Gene.Symbol
probeID <- IDs$Probe.Set.ID

#For matching probes if from different platforms
#  (Section will take ~5-10 minutes to run)
# source("collapseRows_NEW.R") # ONLY uncomment this line if you get an error with it commented
# datExprA5g = (collapseRows(datExprB1,genes,probeID))[[1]]
# datExprA6g = (collapseRows(datExprB2,genes,probeID))[[1]]

#Limit analysis to common probes
commonProbesA = intersect (rownames(datExprA1),rownames(datExprA2))
commonProbesA = intersect (commonProbesA, rownames(datExprA3))
commonProbesA = intersect (commonProbesA, rownames(datExprA4))
commonProbesA = intersect (commonProbesA, rownames(datExprA7))

datExprA1p = datExprA1[commonProbesA,]
datExprA2p = datExprA2[commonProbesA,]
datExprA3p = datExprA3[commonProbesA,]
datExprA4p = datExprA4[commonProbesA,]
datExprA7p = datExprA7[commonProbesA,]

#Convert probes into gene symbols
keepGenesDups = (collapseRows(datExprA1p,genes,probeID))[[2]]
datExprA1g    = datExprA1p[keepGenesDups[,2],]
datExprA2g    = datExprA2p[keepGenesDups[,2],]
datExprA3g    = datExprA3p[keepGenesDups[,2],]
datExprA4g    = datExprA4p[keepGenesDups[,2],]
datExprA7g    = datExprA7p[keepGenesDups[,2],]

rownames(datExprA1g)<-rownames(datExprA2g)<-rownames(datExprA3g)<-rownames(datExprA4g)<-rownames(datExprA7g)<-keepGenesDups[,1]

#Find common genes to all
commonGenesA = intersect (rownames(datExprA1g),rownames(datExprA2g))
commonGenesA = intersect (commonGenesA, rownames(datExprA3g))
commonGenesA = intersect (commonGenesA, rownames(datExprA4g))
commonGenesA = intersect (commonGenesA, rownames(datExprA5g))
commonGenesA = intersect (commonGenesA, rownames(datExprA6g))
commonGenesA = intersect (commonGenesA, rownames(datExprA7g))

#Restrict data set to common genes
datExprA1g    = datExprA1g[commonGenesA,]
datExprA2g    = datExprA2g[commonGenesA,]
datExprA3g    = datExprA3g[commonGenesA,]
datExprA4g    = datExprA4g[commonGenesA,]
datExprA5g    = datExprA5g[commonGenesA,]
datExprA6g    = datExprA6g[commonGenesA,]
datExprA7g    = datExprA7g[commonGenesA,]

# commonGenesB = intersect (rownames(datExprB1g),rownames(datExprB2g))
# datExprB1g = datExprB1g[commonGenesB,]
# datExprB2g = datExprB2g[commonGenesB,]

# Exp <- t(datExprA4g)

# ###Choosing soft threshold
# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(Exp, powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.50,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# #A1 = 7
# #A2 = 9
# #A3 = 14
# #A4 = 9
# #A5 = 9
# #A6 = 9


setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/")
#Correlating general network properties
softPower = 10 # (Read WGCNA tutorial to learn how to pick your power)
rankExprA1= rank(rowMeans(datExprA1g))
rankExprA2= rank(rowMeans(datExprA2g))
rankExprA3= rank(rowMeans(datExprA3g))
rankExprA4= rank(rowMeans(datExprA4g))
rankExprA5= rank(rowMeans(datExprA5g))
rankExprA6= rank(rowMeans(datExprA6g))
rankExprA7= rank(rowMeans(datExprA7g))

random4000= sample(commonGenesA,2800)
rankConnA1= rank(softConnectivity(t(datExprA1g[random4000,]),type="signed",power=softPower, minNSamples = 3))
rankConnA2= rank(softConnectivity(t(datExprA2g[random4000,]),type="signed",power=softPower, minNSamples = 3))
rankConnA3= rank(softConnectivity(t(datExprA3g[random4000,]),type="signed",power=softPower, minNSamples = 3))
rankConnA4= rank(softConnectivity(t(datExprA4g[random4000,]),type="signed",power=softPower, minNSamples = 3))
rankConnA5= rank(softConnectivity(t(datExprA5g[random4000,]),type="signed",power=softPower, minNSamples = 3))
rankConnA6= rank(softConnectivity(t(datExprA6g[random4000,]),type="signed",power=softPower, minNSamples = 3))
rankConnA7= rank(softConnectivity(t(datExprA7g[random4000,]),type="signed",power=softPower, minNSamples = 3))

# rankExprB1= rank(rowMeans(datExprB1g))
# rankExprB2= rank(rowMeans(datExprB2g))
# random5000= sample(commonGenesB,5000)
# rankConnB1= rank(softConnectivity(t(datExprB1g[random5000,]),type="signed",power=softPower, minNSamples = 3))
# rankConnB2= rank(softConnectivity(t(datExprB2g[random5000,]),type="signed",power=softPower, minNSamples = 3))

pdf("generalNetworkProperties.pdf", height=10, width=9)
par(mfrow=c(2,2))
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (C9orf72)", 
                   ylab="Ranked Expression (CHMP2B)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA1,rankConnA2, xlab="Ranked Connectivity (C9orf72)", 
                   ylab="Ranked Connectivity (CHMP2B)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA1,rankExprA3, xlab="Ranked Expression (C9orf72)", 
                   ylab="Ranked Expression (sALS)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA1,rankConnA3, xlab="Ranked Connectivity (C9orf72)", 
                   ylab="Ranked Connectivity (sALS)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA1,rankExprA4, xlab="Ranked Expression (C9orf72)", 
                   ylab="Ranked Expression (VCP)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA1,rankConnA4, xlab="Ranked Connectivity (C9orf72)", 
                   ylab="Ranked Connectivity (VCP)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA1,rankExprA5, xlab="Ranked Expression (C9orf72)", 
                   ylab="Ranked Expression (PET)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA1,rankConnA5, xlab="Ranked Connectivity (C9orf72)", 
                   ylab="Ranked Connectivity (PET)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA1,rankExprA6, xlab="Ranked Expression (C9orf72)", 
                   ylab="Ranked Expression (RAV)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA1,rankConnA6, xlab="Ranked Connectivity (C9orf72)", 
                   ylab="Ranked Connectivity (RAV)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA1,rankExprA7, xlab="Ranked Expression (C9orf72)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA1,rankConnA7, xlab="Ranked Connectivity (C9orf72)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA2,rankExprA3, xlab="Ranked Expression (CHMP2B)", 
                   ylab="Ranked Expression (sALS)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA2,rankConnA3, xlab="Ranked Connectivity (CHMP2B)", 
                   ylab="Ranked Connectivity (sALS)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA2,rankExprA4, xlab="Ranked Expression (CHMP2B)", 
                   ylab="Ranked Expression (VCP)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA2,rankConnA4, xlab="Ranked Connectivity (CHMP2B)", 
                   ylab="Ranked Connectivity (VCP)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA2,rankExprA5, xlab="Ranked Expression (CHMP2B)", 
                   ylab="Ranked Expression (PET)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA2,rankConnA5, xlab="Ranked Connectivity (CHMP2B)", 
                   ylab="Ranked Connectivity (PET)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA2,rankExprA6, xlab="Ranked Expression (CHMP2B)", 
                   ylab="Ranked Expression (RAV)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA2,rankConnA6, xlab="Ranked Connectivity (CHMP2B)", 
                   ylab="Ranked Connectivity (RAV)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA2,rankExprA7, xlab="Ranked Expression (CHMP2B)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA2,rankConnA7, xlab="Ranked Connectivity (CHMP2B)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA3,rankExprA4, xlab="Ranked Expression (sALS)", 
                   ylab="Ranked Expression (VCP)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA3,rankConnA4, xlab="Ranked Connectivity (sALS)", 
                   ylab="Ranked Connectivity (VCP)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA3,rankExprA5, xlab="Ranked Expression (sALS)", 
                   ylab="Ranked Expression (PET)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA3,rankConnA5, xlab="Ranked Connectivity (CHMP2B)", 
                   ylab="Ranked Connectivity (PET)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA3,rankExprA6, xlab="Ranked Expression (sALS)", 
                   ylab="Ranked Expression (RAV)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA3,rankConnA6, xlab="Ranked Connectivity (sALS)", 
                   ylab="Ranked Connectivity (RAV)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA3,rankExprA7, xlab="Ranked Expression (sALS)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA3,rankConnA7, xlab="Ranked Connectivity (sALS)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA4,rankExprA5, xlab="Ranked Expression (VCP)", 
                   ylab="Ranked Expression (PET)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA4,rankConnA5, xlab="Ranked Connectivity (VCP)", 
                   ylab="Ranked Connectivity (PET)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA4,rankExprA6, xlab="Ranked Expression (VCP)", 
                   ylab="Ranked Expression (RAV)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA4,rankConnA6, xlab="Ranked Connectivity (VCP)", 
                   ylab="Ranked Connectivity (RAV)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA4,rankExprA7, xlab="Ranked Expression (VCP)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA4,rankConnA7, xlab="Ranked Connectivity (VCP)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA5,rankExprA6, xlab="Ranked Expression (PET)", 
                   ylab="Ranked Expression (RAV)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA5,rankConnA6, xlab="Ranked Connectivity (PET)", 
                   ylab="Ranked Connectivity (RAV)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA5,rankExprA7, xlab="Ranked Expression (PET)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA5,rankConnA7, xlab="Ranked Connectivity (PET)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprA6,rankExprA7, xlab="Ranked Expression (RAV)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnA6,rankConnA7, xlab="Ranked Connectivity (RAV)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

dev.off()

# verboseScatterplot(rankExprB1,rankExprB2, xlab="Ranked Expression (B1)", 
#                    ylab="Ranked Expression (B2)")
# verboseScatterplot(rankConnB1,rankConnB2, xlab="Ranked Connectivity (B1)", 
#                    ylab="Ranked Connectivity (B2)")



#calculate all of the necessary values to run WGCNA
#(this will take around 10 minutes)
adjacencyA1 = adjacency(t(datExprA1g),power=7,type="signed");
diag(adjacencyA1)=0
dissTOMA1   = 1-TOMsimilarity(adjacencyA1, TOMType="signed")
geneTreeA1  = flashClust(as.dist(dissTOMA1), method="average")

adjacencyA2 = adjacency(t(datExprA2g),power=9,type="signed");
diag(adjacencyA2)=0
dissTOMA2   = 1-TOMsimilarity(adjacencyA2, TOMType="signed")
geneTreeA2  = flashClust(as.dist(dissTOMA2), method="average")

adjacencyA3 = adjacency(t(datExprA3g),power=14,type="signed");
diag(adjacencyA3)=0
dissTOMA3   = 1-TOMsimilarity(adjacencyA3, TOMType="signed")
geneTreeA3  = flashClust(as.dist(dissTOMA3), method="average")

adjacencyA4 = adjacency(t(datExprA4g),power=9,type="signed");
diag(adjacencyA4)=0
dissTOMA4   = 1-TOMsimilarity(adjacencyA4, TOMType="signed")
geneTreeA4  = flashClust(as.dist(dissTOMA4), method="average")

adjacencyA7 = adjacency(t(datExprA7g),power=9,type="signed");
diag(adjacencyA7)=0
dissTOMA7   = 1-TOMsimilarity(adjacencyA7, TOMType="signed")
geneTreeA7  = flashClust(as.dist(dissTOMA7), method="average")



#For RNA seq data, save as a csv file, open in excel and change to "number". Save as txt file and reload

adjacencyA5 = adjacency(t(datExprA5g),power=9,type="signed");
diag(adjacencyA5)=0
dissTOMA5   = 1-TOMsimilarity(adjacencyA5, TOMType="signed")
geneTreeA5  = flashClust(as.dist(dissTOMA5), method="average")

adjacencyA6 = adjacency(t(datExprA6g),power=9,type="signed");
diag(adjacencyA6)=0
dissTOMA6   = 1-TOMsimilarity(adjacencyA6, TOMType="signed")
geneTreeA6  = flashClust(as.dist(dissTOMA6), method="average")

# adjacencyB1 = adjacency(t(datExprB1g),power=softPower,type="signed");
# diag(adjacencyB1)=0
# dissTOMB1   = 1-TOMsimilarity(adjacencyB1, TOMType="signed")
# geneTreeB1  = flashClust(as.dist(dissTOMB1), method="average")
# 
# adjacencyB2 = adjacency(t(datExprB2g),power=softPower,type="signed");
# diag(adjacencyB2)=0
# dissTOMB2   = 1-TOMsimilarity(adjacencyB2, TOMType="signed")
# geneTreeB2  = flashClust(as.dist(dissTOMB2), method="average")

# save.image("tutorial.RData")  #  (Section will take ~5-15 minutes to run)

pdf("dendrogram.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(geneTreeA1,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (C9orf72)", labels=FALSE,hang=0.04);
plot(geneTreeA2,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (CHMP2B)", labels=FALSE,hang=0.04); 
plot(geneTreeA3,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (sALS)", labels=FALSE,hang=0.04);
plot(geneTreeA4,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (VCP)", labels=FALSE,hang=0.04); 
plot(geneTreeA5,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (PET)", labels=FALSE,hang=0.04);
plot(geneTreeA6,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (RAV)", labels=FALSE,hang=0.04);
plot(geneTreeA7,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (FTLD)", labels=FALSE,hang=0.04);
dev.off()

# plot(geneTreeB1,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (B1)", labels=FALSE,hang=0.04);
# plot(geneTreeB2,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (B2)", labels=FALSE,hang=0.04); 


#determine modules based on control data set

mColorh=NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTreeA1, pamStage=FALSE,
                      minClusterSize = (30-3*ds), cutHeight = 0.5, 
                      deepSplit = ds, distM = dissTOMA1)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}
pdf("Module_choices5.pdf", height=10,width=25); 
plotDendroAndColors(geneTreeA1, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);
dev.off()

modulesA1 =  mColorh[,3] #choose based on deepslit values in plot

#calculate the principle components for visualizations 
PCs1A    = moduleEigengenes(t(datExprA1g),  colors=modulesA1) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesA1))
pdf("ModuleEigengeneVisualizations.pdf",height=6,width=8)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)

plot(pcTree1A, xlab="",ylab="",main="",sub="")

plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreeA1$order
plotMat(scale(log(datExprA1g[ordergenes,])) , rlabels= modulesA1[ordergenes], clabels= colnames(datExprA1g), rcols=modulesA1[ordergenes])

for (which.module in names(table(modulesA1))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
} 

dev.off()

##Qualitatively and quantitatively measure network preservation at the module level##

#assess how well modules in network 1 are preserved in network 2
pdf("Final_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTreeA1, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (C9orf72)") 
plotDendroAndColors(geneTreeA2, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (CHMP2B)") 
plotDendroAndColors(geneTreeA3, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (sALS)") 
plotDendroAndColors(geneTreeA4, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (VCP)") 
plotDendroAndColors(geneTreeA5, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (PET)") 
plotDendroAndColors(geneTreeA6, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (RAV)") 

dev.off()


# plotDendroAndColors(geneTreeB1, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
#                     guideHang=0.05, main="Gene dendrogram and module colors (A1)") 
# plotDendroAndColors(geneTreeB2, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
#                     guideHang=0.05, main="Gene dendrogram and module colors (A2)") 




#assess how well a module in one study is preserved in another study
# (This step will take ~10 minutes)
multiExpr  = list(A1=list(data=t(datExprA1g)),A3=list(data=t(datExprA3g)),A4=list(data=t(datExprA4g)),
                  A5=list(data=t(datExprA5g)),A6=list(data=t(datExprA6g)))
multiColor = list(A1 = modulesA1)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400)
stats1 = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A3
stats2 = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A4
stats3 = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A5
stats4 = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A6

stats1[order(-stats1[,2]),c(1:2)]
stats2[order(-stats2[,2]),c(1:2)]
stats3[order(-stats3[,2]),c(1:2)]
stats4[order(-stats4[,2]),c(1:2)]

statslist <- stats[order(-stats[,2]),c(1:2)]

write.csv(stats4, file = "stats4list.csv")

#get the kME values
geneModuleMembership1 = signedKME(t(datExprA1g), ME_1A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExprA1g)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep="");

Gene       = rownames(datExprA1g)
kMEtable1  = cbind(Gene,Gene,modulesA1)
for (i in 1:length(colorsA1))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

write.csv(kMEtable1,"kMEtable1.csv",row.names=FALSE)

# First calculate MEs for A2, since we haven't done that yet
PCs2A = moduleEigengenes(t(datExprA2g),  colors=modulesA1) 
ME_2A = PCs2A$eigengenes

geneModuleMembership2 = signedKME(t(datExprA2g), ME_2A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(datExprA2g)[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsA1,".pval",sep="");

kMEtable2  = cbind(Gene,Gene,modulesA1)
for (i in 1:length(colorsA1))
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)

write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)

#plot the kME values 
pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsA1[c],
                     xlab="kME in A2",ylab="kME in A1")
}; dev.off()

pdf("inModule_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
  inMod = modulesA1== colorsA1[c]
  verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsA1[c],
                     xlab="kME in A2",ylab="kME in A1")
}; dev.off()

#determine which genes are hubs in both networks 
topGenesKME = NULL
for (c in 1:length(colorsA1)){
  kMErank1    = rank(-geneModuleMembership1[,c])
  kMErank2    = rank(-geneModuleMembership2[,c])
  maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsA1
topGenesKME

# ##Comparing networks and annotating modules using programs outside of R##
# 
# #output data from our network for import into VisANT
# source("tutorialFunctions.R")  
# #Dataset 1
# for (co in colorsA1[colorsA1!="grey"])
#   visantPrepOverall(modulesA1, co, t(datExprA1g), rownames(datExprA1g), 500, softPower, TRUE)
# #Dataset 2
# for (co in colorsA1[colorsA1!="grey"])
#   visantPrepOverall(modulesA1, co, t(datExprA2g), rownames(datExprA2g), 500, softPower, TRUE)


# #annotate modules based enrichment for user-defined lists
# 
enrichments = userListEnrichment(Gene, modulesA1, c("TDP43Enrichment.csv","kMEtable1.csv"),
                                 c("enrichmentlist","humanModules"), "enrichment.csv")
write.csv(x = enrichments$pValues, file = "allenrichments.csv")

splitgroups <- split(genemod, genemod[,2])
# 
# #compare how a gene or module relates to phenotype across data sets 
# group <- c("Control","Control","Control","Patient","Patient","Patient","Patient","Patient","Patient","Patient","Patient")
# var <- list(group=="Control", group=="Patient")
# datgroupM = t(apply(t(ME_1A),1,t.test.l))
# colnames(datgroupM)=c("MeanCon","MeanPat","SD_Con","SD_Pat","PvalGroup")
# datgroupM[datgroupM[,5]<0.05,]
# 
# group2 <- c("Control","Control","Control","Control","Control","Control","Patient","Patient","Patient")
# var <- list(group2=="Control", group2=="Patient")
# datgroup2 = t(apply(t(datExprA2g),1,t.test.l))
# colnames(datgroup2)=c("MeanCon","MeanPat","SD_Con","SD_Pat","PvalGroup")
# datgroupM2 = t(apply(t(ME_2A),1,t.test.l))
# colnames(datgroupM2)=c("MeanCon","MeanPat","SD_Con","SD_Pat","PvalGroup")
# 
# datgroupM2[datgroupM2[,5]<0.05,]
