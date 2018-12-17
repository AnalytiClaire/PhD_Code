#### DEG PPI Correlation ####

#Read in network nodes
DEG_PPI <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/OLDPPI/DEG_PPI_Genes_nofib.txt")

#Extract PPI network genes from each dataset
C9 <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/C9uniquegene_samples.csv", row.names = 1)
C9 <- C9[,4:11]
C9 <- subset(C9, rownames(C9) %in% DEG_PPI)

VCP <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/VCPuniquegene_samples.csv", row.names = 1)
VCP <- VCP[,4:10]
VCP <- subset(VCP, rownames(VCP) %in% DEG_PPI)

FTLD <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/FTLDuniquegene_samples.csv", row.names = 1)
FTLD <- FTLD[,9:24]
FTLD <- subset(FTLD, rownames(FTLD) %in% DEG_PPI)

sALS <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/sALSuniquegene_samples.csv", row.names = 1)
sALS <- sALS[,4:10]
sALS <- subset(sALS, rownames(sALS) %in% DEG_PPI)

PET <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/PET_results_keepfiltering.csv")
rownames(PET) <- PET$hgnc_symbol
PET <- PET[,19:35]
PET <- subset(PET, rownames(PET) %in% DEG_PPI)

RAV <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/RAV_results_keepfiltering.csv")
rownames(RAV) <- RAV$hgnc_symbol
RAV <- RAV[,18:30]
RAV <- subset(RAV, rownames(RAV) %in% DEG_PPI)

#Find the gene names that all datasets have in common
DEG_com <- Reduce(intersect, list(rownames(C9), rownames(VCP), rownames(FTLD),
                                  rownames(sALS), rownames(PET), rownames(RAV)))

#Subset each dataset with these common names so they are all the same size
C9 <- subset(C9, rownames(C9) %in% DEG_com)
VCP <- subset(VCP, rownames(VCP) %in% DEG_com)
FTLD <- subset(FTLD, rownames(FTLD) %in% DEG_com)
sALS <- subset(sALS, rownames(sALS) %in% DEG_com)
PET <- subset(PET, rownames(PET) %in% DEG_com)
RAV <- subset(RAV, rownames(RAV) %in% DEG_com)

#### Cor.test Method ####

library(tictoc)
library(gdata)

##For loop for generating regression values and p values
CorExprMat <- t(C9[1:100,])
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
write.csv(reg.mat, "VCP_Rho_PPI_coexpression.csv")


#### filter correlations ####
setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/")

#Split names to two columns
C9_cor <- read.csv("C9_PPI_coexpression.csv")
C9_cor$Gene1 <- as.character(lapply(strsplit(as.character(C9_cor$X), "\\:"), "[", 2))
C9_cor$Gene2 <- as.character(lapply(strsplit(as.character(C9_cor$X), "\\:"), "[", 1))
C9_cor <- C9_cor[,c(5,6,1,2,3,4)]

sALS_cor <- read.csv("sALS_PPI_coexpression.csv")
sALS_cor$Gene1 <- as.character(lapply(strsplit(as.character(sALS_cor$X), "\\:"), "[", 2))
sALS_cor$Gene2 <- as.character(lapply(strsplit(as.character(sALS_cor$X), "\\:"), "[", 1))
sALS_cor <- sALS_cor[,c(5,6,1,2,3,4)]

FTLD_cor <- read.csv("FTLD_PPI_coexpression.csv")
FTLD_cor$Gene1 <- as.character(lapply(strsplit(as.character(FTLD_cor$X), "\\:"), "[", 2))
FTLD_cor$Gene2 <- as.character(lapply(strsplit(as.character(FTLD_cor$X), "\\:"), "[", 1))
FTLD_cor <- FTLD_cor[,c(5,6,1,2,3,4)]

VCP_cor <- read.csv("VCP_PPI_coexpression.csv")
VCP_cor$Gene1 <- as.character(lapply(strsplit(as.character(VCP_cor$X), "\\:"), "[", 2))
VCP_cor$Gene2 <- as.character(lapply(strsplit(as.character(VCP_cor$X), "\\:"), "[", 1))
VCP_cor <- VCP_cor[,c(5,6,2,3,4)]

PET_cor <- read.csv("PET_PPI_coexpression.csv")
PET_cor$Gene1 <- as.character(lapply(strsplit(as.character(PET_cor$X), "\\:"), "[", 2))
PET_cor$Gene2 <- as.character(lapply(strsplit(as.character(PET_cor$X), "\\:"), "[", 1))
PET_cor <- PET_cor[,c(5,6,2,3,4)]

RAV_cor <- read.csv("RAV_PPI_coexpression.csv")
RAV_cor$Gene1 <- as.character(lapply(strsplit(as.character(RAV_cor$X), "\\:"), "[", 2))
RAV_cor$Gene2 <- as.character(lapply(strsplit(as.character(RAV_cor$X), "\\:"), "[", 1))
RAV_cor <- RAV_cor[,c(5,6,2,3,4)]

hist(C9_cor$reg.mat)
hist(sALS_cor$reg.mat)
hist(VCP_cor$reg.mat)
hist(FTLD_cor$reg.mat)
hist(PET_cor$reg.mat)
hist(RAV_cor$reg.mat)


### Filter by r value
C9_cor.5 <- C9_cor[C9_cor$reg.mat > 0.5 | C9_cor$reg.mat < -0.5,]
# sALS_cor.5 <- sALS_cor[sALS_cor$reg.mat > 0.5 | sALS_cor$reg.mat < -0.5,]
FTLD_cor.5 <- FTLD_cor[FTLD_cor$reg.mat > 0.5 | FTLD_cor$reg.mat < -0.5,]
VCP_cor.5 <- VCP_cor[VCP_cor$reg.mat > 0.5 | VCP_cor$reg.mat < -0.5,]
PET_cor.5 <- PET_cor[PET_cor$reg.mat > 0.5 | PET_cor$reg.mat < -0.5,]
RAV_cor.5 <- RAV_cor[RAV_cor$reg.mat > 0.5 | RAV_cor$reg.mat < -0.5,]


Commonedge <- Reduce(intersect, list(
  FTLD_cor.5$X,
  VCP_cor.5$X,
  PET_cor.5$X, 
  RAV_cor.5$
))









setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/")
write.csv(C9_cor.5, "C9_cor.5_cytoscape_combo.csv", row.names = F, quote = F)
write.csv(sALS_cor.5, "sALS_cor.5_cytoscape_combo.csv", row.names = F, quote = F)
write.csv(FTLD_cor.5, "FTLD_cor.5_cytoscape_combo.csv", row.names = F, quote = F)
write.csv(VCP_cor.5, "VCP_cor.5_cytoscape_combo.csv", row.names = F, quote = F)
write.csv(PET_cor.5, "PET_cor.5_cytoscape_combo.csv", row.names = F, quote = F)
write.csv(RAV_cor.5, "RAV_cor.5_cytoscape_combo.csv", row.names = F, quote = F)


#Investigating how cytoscape merges edges
# 
# setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/")
# C9_cor.5 <- read.csv("C9_cor.5_cytoscape.csv")
# sALS_cor.5 <-read.csv("sALS_cor.5_cytoscape.csv")
# FTLD_cor.5 <-read.csv("FTLD_cor.5_cytoscape.csv")
# VCP_cor.5 <-read.csv("VCP_cor.5_cytoscape.csv")
# PET_cor.5 <-read.csv("PET_cor.5_cytoscape.csv")
# RAV_cor.5 <-read.csv("RAV_cor.5_cytoscape.csv")

Group2 <- read.csv("SubnetworkEdgeInfo.csv")

C9merge <- merge(Group2, C9_cor.5, by.x = "X", by.y = "combo")
sALSmerge <- merge(Group2, sALS_cor.5, by.x = "X", by.y = "combo")
FTLDmerge <- merge(Group2, FTLD_cor.5, by.x = "X", by.y = "combo")
VCPmerge <- merge(Group2, VCP_cor.5, by.x = "X", by.y = "combo")
PETmerge <- merge(Group2, PET_cor.5, by.x = "X", by.y = "combo")
RAVmerge <- merge(Group2, RAV_cor.5, by.x = "X", by.y = "combo")

#Save common module as edge file and split gene names into two columns

#Merge into one table. First data is read in, but because the genes have been separated into two columns
#we want them back in one column so that they can be cross referenced correctly. It just so happens that to match
#the cytoscape format, we have to put the gene in the second column first and first column second.

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/")
C9_cor.5 <- read.csv("C9_cor.5_cytoscape.csv")
C9_cor.5$combo <- paste(C9_cor.5$Gene2,":",C9_cor.5$Gene1, sep = "")

sALS_cor.5 <-read.csv("sALS_cor.5_cytoscape.csv")
sALS_cor.5$combo <- paste(sALS_cor.5$Gene2,":",sALS_cor.5$Gene1, sep = "")

FTLD_cor.5 <-read.csv("FTLD_cor.5_cytoscape.csv")
FTLD_cor.5$combo <- paste(FTLD_cor.5$Gene2,":",FTLD_cor.5$Gene1, sep = "")

VCP_cor.5 <-read.csv("VCP_cor.5_cytoscape.csv")
VCP_cor.5$combo <- paste(VCP_cor.5$Gene2,":",VCP_cor.5$Gene1, sep = "")

PET_cor.5 <-read.csv("PET_cor.5_cytoscape.csv")
PET_cor.5$combo <- paste(PET_cor.5$Gene2,":",PET_cor.5$Gene1, sep = "")

RAV_cor.5 <-read.csv("RAV_cor.5_cytoscape.csv")
RAV_cor.5$combo <- paste(RAV_cor.5$Gene2,":",RAV_cor.5$Gene1, sep = "")


#Original module
# Group1 <- read.csv("MergedNetworkdefaultedgeNEW.csv")
#Recalculated module (CORRECT VERSION)
Group2 <- read.csv("SubnetworkEdgeInfo.csv")

C9merge <- merge(Group2, C9_cor.5, by.x = "X", by.y = "combo")
sALSmerge <- merge(Group2, sALS_cor.5, by.x = "X", by.y = "combo")
FTLDmerge <- merge(Group2, FTLD_cor.5, by.x = "X", by.y = "combo")
VCPmerge <- merge(Group2, VCP_cor.5, by.x = "X", by.y = "combo")
PETmerge <- merge(Group2, PET_cor.5, by.x = "X", by.y = "combo")
RAVmerge <- merge(Group2, RAV_cor.5, by.x = "X", by.y = "combo")


Group1Expr <- data.frame(row.names = C9merge$name,
                         C9 = C9merge$reg.mat.y,
                         sALS = sALSmerge$reg.mat.y,
                         FTLD = FTLDmerge$reg.mat.y,
                         VCP = VCPmerge$reg.mat.y,
                         PET = PETmerge$reg.mat.y,
                         RAV = RAVmerge$reg.mat.y)

Group1up <- subset(Group1Expr, (Group1Expr$C9 > 0.5 &
                                  Group1Expr$sALS > 0.5 &
                                  Group1Expr$FTLD > 0.5 &
                                  Group1Expr$VCP > 0.5 &
                                  Group1Expr$PET > 0.5 &
                                  Group1Expr$RAV > 0.5))

Group1down <- subset(Group1Expr, (Group1Expr$C9 < -0.5 &
                                    Group1Expr$sALS < -0.5 &
                                    Group1Expr$FTLD < -0.5 &
                                    Group1Expr$VCP < -0.5 &
                                    Group1Expr$PET < -0.5 &
                                    Group1Expr$RAV < -0.5))

Group1samedir <- rbind(Group1up, Group1down)
Group1samedir$corMean <- rowMeans(Group1samedir, na.rm = FALSE, dims = 1)
write.csv(Group1samedir, "CORRECTEDG1_samedir_mean.csv", quote = F)

#### Overlap with DEGs ####
DEGs <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_upanddown.txt")
ALSOD <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/ALSoDgenes.txt")
Disgenes <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")

overlap <- Reduce(intersect, list(Disgenes, Group1Genes))


#### Repeat for controls ####

#Read in network nodes
Group1Genes <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NuggetGenes.txt")

#Extract PPI network genes from each dataset
C9 <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/C9uniquegene_samples.csv", row.names = 1)
C9 <- C9[,1:3]
C9 <- subset(C9, rownames(C9) %in% Group1Genes)

VCP <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/VCPuniquegene_samples.csv", row.names = 1)
VCP <- VCP[,1:3]
VCP <- subset(VCP, rownames(VCP) %in% Group1Genes)

FTLD <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/FTLDuniquegene_samples.csv", row.names = 1)
FTLD <- FTLD[,1:8]
FTLD <- subset(FTLD, rownames(FTLD) %in% Group1Genes)

sALS <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/sALSuniquegene_samples.csv", row.names = 1)
sALS <- sALS[,1:3]
sALS <- subset(sALS, rownames(sALS) %in% Group1Genes)

PET <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/PET_results_keepfiltering.csv")
rownames(PET) <- PET$hgnc_symbol
PET <- PET[,10:18]
PET <- subset(PET, rownames(PET) %in% Group1Genes)

RAV <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/RAV_results_keepfiltering.csv")
rownames(RAV) <- RAV$hgnc_symbol
RAV <- RAV[,10:17]
RAV <- subset(RAV, rownames(RAV) %in% Group1Genes)


#### Cor.test Method ####

library(tictoc)
library(gdata)

# ##For loop for generating regression values and p values
# CorExprMat <- t(RAV)
# 
# test <- CorExprMat
# 
# reg <- matrix(0, ncol(test), ncol(test))
# p.value <- matrix(0, ncol(test), ncol(test))
# 
# tic()
# for (i in 1:ncol(test)){
#   for (j in 1:ncol(test)){
#     reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
#   }}
# 
# rownames(reg) <- colnames(reg) <- colnames(test)
# toc()
# 
# tic()
# for (i in 1:ncol(test)){
#   for (j in 1:ncol(test)){
#     p.value[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$p.value
#   }}
# 
# rownames(p.value) <- colnames(p.value) <- colnames(test)
# toc()
# 
# 
# ##Only take upper triangle without diagonal (all comparisons are currently doubled)
# ptri <- p.value
# ptri[lower.tri(ptri, diag = TRUE)] <- NA
# 
# #Turn into vector
# p.vec <- unmatrix(ptri)
# #Remove NA values
# p.vec <- na.omit(p.vec)
# #Multiple hypothesis testing correction
# p.adj <- p.adjust(p.vec, method = "fdr", n = length(p.vec))
# 
# #Create results table
# reg.mat <- unmatrix(reg)
# reg.mat <- as.data.frame(reg.mat)
# p.adj <- as.data.frame(p.adj)
# p.mat <- as.data.frame(p.vec)
# 
# pvals <- merge(p.adj, p.mat, by.x = "row.names", by.y = "row.names")
# rownames(pvals)<- pvals$Row.names
# pvals[,1] <- NULL
# results <- merge(pvals, reg.mat, by.x = "row.names", by.y = "row.names")
# rownames(results)<- results$Row.names
# results[,1] <- NULL
# results <- results[order(results$p.vec),]
# 
# setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/")
# write.csv(results, "RAVcon_PPI_coexpression.csv")


#### Split Names ####

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/")
C9_cor <- read.csv("C9con_PPI_coexpression.csv")
sALS_cor <-read.csv("sALScon_PPI_coexpression.csv")
FTLD_cor <-read.csv("FTLDcon_PPI_coexpression.csv")
VCP_cor <-read.csv("VCPcon_PPI_coexpression.csv")
PET_cor <-read.csv("PETcon_PPI_coexpression.csv")
RAV_cor <-read.csv("RAVcon_PPI_coexpression.csv")

C9_cor.5 <- C9_cor[C9_cor$reg.mat > 0.5 | C9_cor$reg.mat < -0.5,]
sALS_cor.5 <- sALS_cor[sALS_cor$reg.mat > 0.5 | sALS_cor$reg.mat < -0.5,]
FTLD_cor.5 <- FTLD_cor[FTLD_cor$reg.mat > 0.5 | FTLD_cor$reg.mat < -0.5,]
VCP_cor.5 <- VCP_cor[VCP_cor$reg.mat > 0.5 | VCP_cor$reg.mat < -0.5,]
PET_cor.5 <- PET_cor[PET_cor$reg.mat > 0.5 | PET_cor$reg.mat < -0.5,]
RAV_cor.5 <- RAV_cor[RAV_cor$reg.mat > 0.5 | RAV_cor$reg.mat < -0.5,]

Commonedge <- Reduce(intersect, list(C9_cor.5$X, sALS_cor.5$X, FTLD_cor.5$X,
                                     VCP_cor.5$X, PET_cor.5$X, RAV_cor.5$X))


#Subset each dataset with these common names so they are all the same size
C9_CE <- subset(C9_cor.5, C9_cor.5$X %in% Commonedge)
C9_CE <- C9_CE[order(C9_CE$X),]
sALS_CE <- subset(sALS_cor.5, sALS_cor.5$X %in% Commonedge)
sALS_CE <- sALS_CE[order(sALS_CE$X),]
FTLD_CE <- subset(FTLD_cor.5, FTLD_cor.5$X %in% Commonedge)
FTLD_CE <- FTLD_CE[order(FTLD_CE$X),]
VCP_CE <- subset(VCP_cor.5, VCP_cor.5$X %in% Commonedge)
VCP_CE <- VCP_CE[order(VCP_CE$X),]
PET_CE <- subset(PET_cor.5, PET_cor.5$X %in% Commonedge)
PET_CE <- PET_CE[order(PET_CE$X),]
RAV_CE <- subset(RAV_cor.5, RAV_cor.5$X %in% Commonedge)
RAV_CE <- RAV_CE[order(RAV_CE$X),]



CommonGroup <- data.frame(row.names = C9_CE$X,
                          C9 = C9_CE$reg.mat,
                          sALS = sALS_CE$reg.mat,
                          FTLD = FTLD_CE$reg.mat,
                          VCP = VCP_CE$reg.mat,
                          PET = PET_CE$reg.mat,
                          RAV = RAV_CE$reg.mat)


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
write.csv(CG_samedir, "RMETHOD_samedir_mean.csv", quote = F)


#### ENRICHMENTS ####
Nodes <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/CORRECTEDG1/CORRECTEDG1node.csv")
Nodes <- Nodes$name
DEGPPI <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes_nofib.txt")
DEGs <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_upanddown.txt")
ALSOD <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/ALSoDgenes.txt")
Disgenes <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")
Taylor <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/Taylor_TDP43.txt")

intersect(Taylor,Nodes)
