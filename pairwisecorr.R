#Selecting DEGS from expression matrix

#Load list of interesting genes
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")
Genelist <- read.csv("DEG_PPI_Genes.txt", header = F)
Genelist <- Genelist$V1

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
salspat <- salsgene[,11:18]
FTLDpat <- FTLDgene[,16:31]
VCPpat <- VCPgene[,11:18]
PETpat <- PETgene[,19:35]
RAVpat <- RAVgene[,18:30]

# #Make gene symbol row names
# rownames(exprs) <- exprs$Gene.Symbol
# exprspat <- exprs[,52:59]

# #Make gene symbol a column
# exprspat <- cbind(exprspat, exprs$Gene.Symbol)
# colnames(exprspat)[length(exprspat)] <- "Gene.Symbol"
# 
# #Merge by interesting gene names with expression to form matrix
# patgene <- merge(Genelist, exprspat, by.x = "Gene", by.y = "Gene.Symbol")
# rownames(patgene) <- patgene$Gene
# patgene[,1] <- NULL
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/Co-expression/")
# write.csv(patgene, file = "C9_DEG_Exprs.csv")



#### Cor.test Method ####

library(tictoc)
library(gdata)

# # #load dataset
# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/WGCNA/C9orf72")
# uniqueresult <- read.csv("C9result.csv")
# rownames(uniqueresult) <- uniqueresult[,1]
# uniqueresult[,1] <- NULL


##For loop for generating regression values and p values
CorExprMat <- t(C9pat)

test <- CorExprMat

reg <- matrix(0, ncol(test), ncol(test))
p.value <- matrix(0, ncol(test), ncol(test))

tic()
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(test)
toc()

tic()
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    p.value[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$p.value
  }}

rownames(p.value) <- colnames(p.value) <- colnames(test)
toc()


##Only take upper triangle without diagonal (all comparisons are currently doubled)
ptri <- p.value
ptri[lower.tri(ptri, diag = TRUE)] <- NA

#Turn into vector
p.vec <- unmatrix(ptri)
#Remove NA values
p.vec <- na.omit(p.vec)
#Multiple hypothesis testing correction
p.adj <- p.adjust(p.vec, method = "fdr", n = length(p.vec))

#Create results table
reg.mat <- unmatrix(reg)
reg.mat <- as.data.frame(reg.mat)
p.adj <- as.data.frame(p.adj)
p.mat <- as.data.frame(p.vec)

pvals <- merge(p.adj, p.mat, by.x = "row.names", by.y = "row.names")
rownames(pvals)<- pvals$Row.names
pvals[,1] <- NULL
results <- merge(pvals, reg.mat, by.x = "row.names", by.y = "row.names")
rownames(results)<- results$Row.names
results[,1] <- NULL
results <- results[order(results$p.vec),]

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/")
write.csv(results, "C9_PPI_coexpression.csv")



#####
###PSYCH METHOD### 
library(psych)
library(tictoc)
library(gdata)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/WGCNA/C9orf72")

#Generate matrix of correlation and p values
genenames <- read.csv(file = "/Users/clairegreen/Desktop/Allconsensus.txt")
genenames <- genenames[,1]

Cor8000 <- uniqueresult[row.names(uniqueresult) %in% genenames,]
tCor8000 <- t(Cor8000)

tic()
cortest <- corr.test(tCor8000, use = "pairwise", method = "spearman", adjust = "fdr") #must use transposed matrix (genes are colnames)
toc()

#Extract R values
cortestoutput <- cortest$r
corRadj <- cortestoutput
corRadj[lower.tri(corRadj, diag = TRUE)] <- NA

#Turn into vector
corRadj <- as.matrix(corRadj)
corRvec <- unmatrix(corRadj)
#Remove NA values
corRvec <- na.omit(corRvec)
corRvec <- as.data.frame(corRvec)


#Extract P values
cortestpadjust <- cortest$p
corPadj <- cortestpadjust
corPadj[lower.tri(corPadj, diag = TRUE)] <- NA

#Turn into vector
corPadj <- as.matrix(corPadj)
corPvec <- unmatrix(corPadj)
#Remove NA values
corPvec <- na.omit(corPvec)
corPvec <- as.data.frame(corPvec)


CorData <- merge(corRvec, corPvec, by.x = "row.names", by.y = "row.names")

#Select significant results
pointone.output <- subset(CorData, CorData$corPvec < 0.1)
write.csv(CorData, file = "allcor.csv")
write.csv(pointone.output, file = "pointone_cor.csv")


