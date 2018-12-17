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


PETExpr <- read.csv(file = "PETEXPRSrankeduniqueresult.csv")
rownames(PETExpr) <- PETExpr$hgnc_symbol
datExprPET <- PETExpr[,9:26]

RAVExpr <- read.csv(file = "RAVEXPRSrankeduniqueresult.csv")
rownames(RAVExpr) <- RAVExpr$hgnc_symbol
datExprRAV <- RAVExpr[,17:28]

FTLDExpr <- read.csv(file = "FTLDrankeduniqueresult.csv")
rownames(FTLDExpr) <- FTLDExpr$Probe.Set.ID
datExprFTLD <- FTLDExpr[,57:72]
annotation <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/GoodData/FTLDID geneIDs.csv")
datExprFTLD <- merge(datExprFTLD, annotation, by.x = "row.names", by.y = "affy_hg_u133_plus_2")
datExprFTLD$ensembl_gene_id <- NULL
datExprFTLD <- datExprFTLD[!duplicated(datExprFTLD[,18]),]
row.names(datExprFTLD) <- datExprFTLD$hgnc_symbol
datExprFTLD$Row.names <- NULL
datExprFTLD$hgnc_symbol <- NULL

#Find common genes to all
commonGenesA = intersect (rownames(datExprPET),rownames(datExprRAV))
commonGenesA = intersect (commonProbesA, rownames(datExprFTLD))

#Restrict data set to common genes
datExprPET    = datExprPET[commonGenesA,]
datExprRAV    = datExprRAV[commonGenesA,]
datExprFTLD    = datExprFTLD[commonGenesA,]


Exp <- t(datExprPET)

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

#FTLD = 5
#RAV = 6
#PET = 6

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/GoodData/")
#Correlating general network properties
softPower = 6 # (Read WGCNA tutorial to learn how to pick your power)
rankExprPET= rank(rowMeans(datExprPET))
rankExprRAV= rank(rowMeans(datExprRAV))
rankExprFTLD= rank(rowMeans(datExprFTLD))


random5000= sample(commonGenesA,7000)
rankConnPET= rank(softConnectivity(t(datExprPET[random5000,]),type="signed",power=softPower))
rankConnRAV= rank(softConnectivity(t(datExprRAV[random5000,]),type="signed",power=softPower))
rankConnFTLD= rank(softConnectivity(t(datExprFTLD[random5000,]),type="signed",power=softPower))

pdf("generalNetworkProperties.pdf", height=10, width=9)
par(mfrow=c(2,2))
verboseScatterplot(rankExprPET,rankExprRAV, xlab="Ranked Expression (PET)", 
                   ylab="Ranked Expression (RAV)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnPET,rankConnRAV, xlab="Ranked Connectivity (PET)", 
                   ylab="Ranked Connectivity (RAV)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprPET,rankExprFTLD, xlab="Ranked Expression (PET)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnPET,rankConnFTLD, xlab="Ranked Connectivity (PET)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

verboseScatterplot(rankExprRAV,rankExprFTLD, xlab="Ranked Expression (RAV)", 
                   ylab="Ranked Expression (FTLD)", abline = TRUE, abline.color = "red")
verboseScatterplot(rankConnRAV,rankConnFTLD, xlab="Ranked Connectivity (RAV)", 
                   ylab="Ranked Connectivity (FTLD)", abline = TRUE, abline.color = "red")

dev.off()

#For RNA seq data, save as a csv file, open in excel and change to "number". Save as txt file and reload
#May have to give rowname column a name

write.table(datExprPET, file = "datExprPET.txt")
write.table(datExprRAV, file = "datExprRAV.txt")
datExprPET <- read.table("datExprPET.txt", header = TRUE)
rownames(datExprPET) <- datExprPET$Gene
datExprRAV <- read.table("datExprRAV.txt", header = TRUE)


#calculate all of the necessary values to run WGCNA
#(this will take around 10 minutes)
adjacencyPET = adjacency(t(datExprPET),power=7,type="signed");
diag(adjacencyPET)=0
dissTOMPET   = 1-TOMsimilarity(adjacencyPET, TOMType="signed")
geneTreePET  = flashClust(as.dist(dissTOMPET), method="average")

adjacencyRAV = adjacency(t(datExprRAV),power=9,type="signed");
diag(adjacencyRAV)=0
dissTOMRAV   = 1-TOMsimilarity(adjacencyRAV, TOMType="signed")
geneTreeRAV  = flashClust(as.dist(dissTOMRAV), method="average")

adjacencyFTLD = adjacency(t(datExprFTLD),power=14,type="signed");
diag(adjacencyFTLD)=0
dissTOMFTLD   = 1-TOMsimilarity(adjacencyFTLD, TOMType="signed")
geneTreeFTLD  = flashClust(as.dist(dissTOMFTLD), method="average")



