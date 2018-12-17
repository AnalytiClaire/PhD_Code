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
PET_cor <- PET_cor[,c(5,6,1,2,3,4)]

RAV_cor <- read.csv("RAV_PPI_coexpression.csv")
RAV_cor$Gene1 <- as.character(lapply(strsplit(as.character(RAV_cor$X), "\\:"), "[", 2))
RAV_cor$Gene2 <- as.character(lapply(strsplit(as.character(RAV_cor$X), "\\:"), "[", 1))
RAV_cor <- RAV_cor[,c(5,6,1,2,3,4)]

hist(C9_cor$reg.mat, main = "C9orf72", xlab = "Correlation Value")
hist(sALS_cor$reg.mat, main = "sALS", xlab = "Correlation Value")
hist(FTLD_cor$reg.mat, main = "FTLD", xlab = "Correlation Value")
hist(VCP_cor$reg.mat, main = "VCP", xlab = "Correlation Value")
hist(PET_cor$reg.mat, main = "PET", xlab = "Correlation Value")
hist(RAV_cor$reg.mat, main = "RAV", xlab = "Correlation Value")

thresh <- 0.55
### Filter by r value
C9_cor.5 <- C9_cor[C9_cor$reg.mat > thresh | C9_cor$reg.mat < -thresh,]
sALS_cor.5 <- sALS_cor[sALS_cor$reg.mat > thresh | sALS_cor$reg.mat < -thresh,]
FTLD_cor.5 <- FTLD_cor[FTLD_cor$reg.mat > thresh | FTLD_cor$reg.mat < -thresh,]
VCP_cor.5 <- VCP_cor[VCP_cor$reg.mat > thresh | VCP_cor$reg.mat < -thresh,]
PET_cor.5 <- PET_cor[PET_cor$reg.mat > thresh | PET_cor$reg.mat < -thresh,]
RAV_cor.5    <- RAV_cor[RAV_cor$reg.mat > thresh | RAV_cor$reg.mat < -thresh,]


C9_cor.5$X2 <- paste(C9_cor.5$Gene1,":",C9_cor.5$Gene2, sep = "")
sALS_cor.5$X2 <- paste(sALS_cor.5$Gene1,":",sALS_cor.5$Gene2, sep = "")
FTLD_cor.5$X2 <- paste(FTLD_cor.5$Gene1,":",FTLD_cor.5$Gene2, sep = "")
VCP_cor.5$X2 <- paste(VCP_cor.5$Gene1,":",VCP_cor.5$Gene2, sep = "")
PET_cor.5$X2 <- paste(PET_cor.5$Gene1,":",PET_cor.5$Gene2, sep = "")
RAV_cor.5$X2 <- paste(RAV_cor.5$Gene1,":",RAV_cor.5$Gene2, sep = "")



#### FULL MERGE ####
#Merge first two datasets
corresult1 <- merge(C9_cor.5, sALS_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(C9_cor.5, sALS_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sALS")

#Merge with 3rd Dataset
corresult3 <- merge(result, FTLD_cor.5, by.x = "Correlation", by.y = "X")
corresult4 <- merge(result, FTLD_cor.5, by.x = "Correlation", by.y = "X2")

corresult3 <- corresult3[, c(1:5, 10)]
corresult4 <- corresult4[, c(1:5, 11)]

result <- rbind(corresult3, corresult4)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sALS", "FTLD")

#Merge with 4th Dataset
corresult5 <- merge(result, VCP_cor.5, by.x = "Correlation", by.y = "X")
corresult6 <- merge(result, VCP_cor.5, by.x = "Correlation", by.y = "X2")

corresult5 <- corresult5[, c(1:6, 7)]
corresult6 <- corresult6[, c(1:6, 8)]

result <- rbind(corresult5, corresult6)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sALS", "FTLD", "VCP")

#Merge with 5th Dataset
corresult7 <- merge(result, PET_cor.5, by.x = "Correlation", by.y = "X")
corresult8 <- merge(result, PET_cor.5, by.x = "Correlation", by.y = "X2")

corresult7 <- corresult7[, c(1:7, 12)]
corresult8 <- corresult8[, c(1:7, 13)]

result <- rbind(corresult7, corresult8)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sALS", "FTLD", "VCP", "PET")

#Merge with 6th Dataset
corresult9 <- merge(result, RAV_cor.5, by.x = "Correlation", by.y = "X")
corresult10 <- merge(result, RAV_cor.5, by.x = "Correlation", by.y = "X2")

corresult9 <- corresult9[, c(1:8, 13)]
corresult10 <- corresult10[, c(1:8, 14)]

result <- rbind(corresult9, corresult10)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sALS", "FTLD", "VCP", "PET", "RAV")

CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:9]


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "Newmethod_5point55.csv", quote = F, row.names = F)

#### Without sALS ####

#Merge first two datasets
corresult1 <- merge(C9_cor.5, FTLD_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(C9_cor.5, FTLD_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD")

#Merge with 3rd Dataset
corresult3 <- merge(result, VCP_cor.5, by.x = "Correlation", by.y = "X")
corresult4 <- merge(result, VCP_cor.5, by.x = "Correlation", by.y = "X2")

corresult3 <- corresult3[, c(1:5, 6)]
corresult4 <- corresult4[, c(1:5, 7)]

result <- rbind(corresult3, corresult4)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD", "VCP")

#Merge with 4th Dataset
corresult5 <- merge(result, PET_cor.5, by.x = "Correlation", by.y = "X")
corresult6 <- merge(result, PET_cor.5, by.x = "Correlation", by.y = "X2")

corresult5 <- corresult5[, c(1:6, 11)]
corresult6 <- corresult6[, c(1:6, 12)]

result <- rbind(corresult5, corresult6)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD", "VCP", "PET")

#Merge with 5th Dataset
corresult7 <- merge(result, RAV_cor.5, by.x = "Correlation", by.y = "X")
corresult8 <- merge(result, RAV_cor.5, by.x = "Correlation", by.y = "X2")

corresult7 <- corresult7[, c(1:7, 12)]
corresult8 <- corresult8[, c(1:7, 13)]

result <- rbind(corresult7, corresult8)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD", "VCP", "PET", "RAV")



CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:8]


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "Newmethod_5point55.csv", quote = F, row.names = F)



### Repeat for Controls ####
#Read in network nodes
Group1Genes <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/NuggetGenes.txt")

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

# sALS <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/sALSuniquegene_samples.csv", row.names = 1)
# sALS <- sALS[,1:3]
# sALS <- subset(sALS, rownames(sALS) %in% Group1Genes)

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

##For loop for generating regression values and p values
CorExprMat <- t(C9)

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

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/")
write.csv(results, "C9con_PPI_coexpression.csv")


#### filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/")
#Split names to two columns
C9_cor <- read.csv("C9con_PPI_coexpression.csv")
C9_cor$Gene1 <- as.character(lapply(strsplit(as.character(C9_cor$X), "\\:"), "[", 2))
C9_cor$Gene2 <- as.character(lapply(strsplit(as.character(C9_cor$X), "\\:"), "[", 1))
C9_cor <- C9_cor[,c(5,6,1,2,3,4)]

FTLD_cor <- read.csv("FTLDcon_PPI_coexpression.csv")
FTLD_cor$Gene1 <- as.character(lapply(strsplit(as.character(FTLD_cor$X), "\\:"), "[", 2))
FTLD_cor$Gene2 <- as.character(lapply(strsplit(as.character(FTLD_cor$X), "\\:"), "[", 1))
FTLD_cor <- FTLD_cor[,c(5,6,1,2,3,4)]

VCP_cor <- read.csv("VCPcon_PPI_coexpression.csv")
VCP_cor$Gene1 <- as.character(lapply(strsplit(as.character(VCP_cor$X), "\\:"), "[", 2))
VCP_cor$Gene2 <- as.character(lapply(strsplit(as.character(VCP_cor$X), "\\:"), "[", 1))
VCP_cor <- VCP_cor[,c(5,6,1,2,3,4)]

PET_cor <- read.csv("PETcon_PPI_coexpression.csv")
PET_cor$Gene1 <- as.character(lapply(strsplit(as.character(PET_cor$X), "\\:"), "[", 2))
PET_cor$Gene2 <- as.character(lapply(strsplit(as.character(PET_cor$X), "\\:"), "[", 1))
PET_cor <- PET_cor[,c(5,6,1,2,3,4)]

RAV_cor <- read.csv("RAVcon_PPI_coexpression.csv")
RAV_cor$Gene1 <- as.character(lapply(strsplit(as.character(RAV_cor$X), "\\:"), "[", 2))
RAV_cor$Gene2 <- as.character(lapply(strsplit(as.character(RAV_cor$X), "\\:"), "[", 1))
RAV_cor <- RAV_cor[,c(5,6,1,2,3,4)]

thresh <- 0.55
### Filter by r value
C9_cor.5 <- C9_cor[C9_cor$reg.mat > thresh | C9_cor$reg.mat < -thresh,]
FTLD_cor.5 <- FTLD_cor[FTLD_cor$reg.mat > thresh | FTLD_cor$reg.mat < -thresh,]
VCP_cor.5 <- VCP_cor[VCP_cor$reg.mat > thresh | VCP_cor$reg.mat < -thresh,]
PET_cor.5 <- PET_cor[PET_cor$reg.mat > thresh | PET_cor$reg.mat < -thresh,]
RAV_cor.5    <- RAV_cor[RAV_cor$reg.mat > thresh | RAV_cor$reg.mat < -thresh,]


C9_cor.5$X2 <- paste(C9_cor.5$Gene1,":",C9_cor.5$Gene2, sep = "")
FTLD_cor.5$X2 <- paste(FTLD_cor.5$Gene1,":",FTLD_cor.5$Gene2, sep = "")
VCP_cor.5$X2 <- paste(VCP_cor.5$Gene1,":",VCP_cor.5$Gene2, sep = "")
PET_cor.5$X2 <- paste(PET_cor.5$Gene1,":",PET_cor.5$Gene2, sep = "")
RAV_cor.5$X2 <- paste(RAV_cor.5$Gene1,":",RAV_cor.5$Gene2, sep = "")

#Merge first two datasets
corresult1 <- merge(C9_cor.5, FTLD_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(C9_cor.5, FTLD_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD")

#Merge with 3rd Dataset
corresult3 <- merge(result, VCP_cor.5, by.x = "Correlation", by.y = "X")
corresult4 <- merge(result, VCP_cor.5, by.x = "Correlation", by.y = "X2")

corresult3 <- corresult3[, c(1:5, 10)]
corresult4 <- corresult4[, c(1:5, 11)]

result <- rbind(corresult3, corresult4)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD", "VCP")

#Merge with 4th Dataset
corresult5 <- merge(result, PET_cor.5, by.x = "Correlation", by.y = "X")
corresult6 <- merge(result, PET_cor.5, by.x = "Correlation", by.y = "X2")

corresult5 <- corresult5[, c(1:6, 11)]
corresult6 <- corresult6[, c(1:6, 12)]

result <- rbind(corresult5, corresult6)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD", "VCP", "PET")

#Merge with 5th Dataset
corresult7 <- merge(result, RAV_cor.5, by.x = "Correlation", by.y = "X")
corresult8 <- merge(result, RAV_cor.5, by.x = "Correlation", by.y = "X2")

corresult7 <- corresult7[, c(1:7, 12)]
corresult8 <- corresult8[, c(1:7, 13)]

result <- rbind(corresult7, corresult8)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "FTLD", "VCP", "PET", "RAV")



CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:8]


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "Con_5point55.csv", quote = F, row.names = F)


##### Repeat for SOD1 and FUS datasets ######

#Read in network nodes
Group1Genes <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/TDP_fixed_NuggetGenes.txt")

#Extract PPI network genes from each dataset
SOD1 <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/non-TDP/SOD1rankeduniqueresult.csv", row.names = 1)
rownames(SOD1) <- SOD1$Gene.Symbol
SOD1 <- SOD1[,27:29]
SOD1 <- subset(SOD1, rownames(SOD1) %in% Group1Genes)

FUS <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/non-TDP/FUSrankeduniqueresult.csv", row.names = 1)
rownames(FUS) <- FUS$Gene.Symbol
FUS <- FUS[,23:25]
FUS <- subset(FUS, rownames(FUS) %in% Group1Genes)


#### Cor.test Method ####

library(tictoc)
library(gdata)

##For loop for generating regression values and p values
CorExprMat <- t(FUS)

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

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/")
write.csv(results, "FUS_PPI_coexpression.csv")


#### filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/")
#Split names to two columns
SOD1_cor <- read.csv("SOD1_PPI_coexpression.csv")
SOD1_cor$Gene1 <- as.character(lapply(strsplit(as.character(SOD1_cor$X), "\\:"), "[", 2))
SOD1_cor$Gene2 <- as.character(lapply(strsplit(as.character(SOD1_cor$X), "\\:"), "[", 1))
SOD1_cor <- SOD1_cor[,c(5,6,1,2,3,4)]

FUS_cor <- read.csv("FUS_PPI_coexpression.csv")
FUS_cor$Gene1 <- as.character(lapply(strsplit(as.character(FUS_cor$X), "\\:"), "[", 2))
FUS_cor$Gene2 <- as.character(lapply(strsplit(as.character(FUS_cor$X), "\\:"), "[", 1))
FUS_cor <- FUS_cor[,c(5,6,1,2,3,4)]

hist(SOD1_cor$reg.mat)
hist(FUS_cor$reg.mat)

thresh <- 0.55
### Filter by r value
SOD1_cor.5 <- SOD1_cor[SOD1_cor$reg.mat > thresh | SOD1_cor$reg.mat < -thresh,]
FUS_cor.5 <- FUS_cor[FUS_cor$reg.mat > thresh | FUS_cor$reg.mat < -thresh,]



SOD1_cor.5$X2 <- paste(SOD1_cor.5$Gene1,":",SOD1_cor.5$Gene2, sep = "")
FUS_cor.5$X2 <- paste(FUS_cor.5$Gene1,":",FUS_cor.5$Gene2, sep = "")


#Merge first two datasets
corresult1 <- merge(SOD1_cor.5, FUS_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(SOD1_cor.5, FUS_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "SOD1", "FUS")



CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:5]


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "SOD1FUS_5point55.csv", quote = F, row.names = F)
