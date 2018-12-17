#### Parkinson's DEG PPI Correlation ####

#Read in network nodes

#### Blood Only ####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/PPI/bloodonly_genes.txt")

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
#Extract PPI network genes from each dataset
LEW <- read.csv("LEWfilteredresult.csv")
rownames(LEW) <- LEW$Gene.Symbol
LEW <- LEW[,28:41]
LEW <- subset(LEW, rownames(LEW) %in% DEG_PPI)

MID3 <- read.csv("MID3filteredresult.csv")
rownames(MID3) <- MID3$Gene.Symbol
MID3 <- MID3[,36:49]
MID3 <- subset(MID3, rownames(MID3) %in% DEG_PPI)

MID4 <- read.csv("MID4filteredresult.csv")
rownames(MID4) <- MID4$Gene.Symbol
MID4 <- MID4[,41:55]
MID4 <- subset(MID4, rownames(MID4) %in% DEG_PPI)

MOR.FC <- read.csv("MOR.FCfilteredresult.csv")
rownames(MOR.FC) <- MOR.FC$Gene.Symbol
MOR.FC <- MOR.FC[,24:28]
MOR.FC <- subset(MOR.FC, rownames(MOR.FC) %in% DEG_PPI)

DIJ <- read.csv("DIJfilteredresult.csv")
rownames(DIJ) <- DIJ$Gene.Symbol
DIJ <- DIJ[,29:43]
DIJ <- subset(DIJ, rownames(DIJ) %in% DEG_PPI)

FFR <- read.csv("FFRfilteredresult.csv")
rownames(FFR) <- FFR$Gene.Symbol
FFR <- FFR[,29:44]
FFR <- subset(FFR, rownames(FFR) %in% DEG_PPI)

MID1 <- read.csv("MID1filteredresult.csv")
rownames(MID1) <- MID1$Gene.Symbol
MID1 <- MID1[,28:37]
MID1 <- subset(MID1, rownames(MID1) %in% DEG_PPI)

MID2 <- read.csv("MID2filteredresult.csv")
rownames(MID2) <- MID2$Gene.Symbol
MID2 <- MID2[,39:49]
MID2 <- subset(MID2, rownames(MID2) %in% DEG_PPI)

MOR.SN <- read.csv("MOR.SNfilteredresult.csv")
rownames(MOR.SN) <- MOR.SN$Gene.Symbol
MOR.SN <- MOR.SN[,36:59]
MOR.SN <- subset(MOR.SN, rownames(MOR.SN) %in% DEG_PPI)

DUM <- read.csv("DUM_UniqueGene_DESeq2.csv")
rownames(DUM) <- DUM$hgnc_symbol
DUM <- DUM[,53:81]
DUM <- subset(DUM, rownames(DUM) %in% DEG_PPI)

CON <- read.csv("CONfilteredresult.csv")
rownames(CON) <- CON$Gene.Symbol
CON <- CON[,29:34]
CON <- subset(CON, rownames(CON) %in% DEG_PPI)

BOT <- read.csv("BOTrankeduniqueresult.csv")
rownames(BOT) <- BOT$Gene.Symbol
BOT <- BOT[,13:14]
BOT <- subset(BOT, rownames(BOT) %in% DEG_PPI)

BOT2 <- read.csv("BOT2rankeduniqueresult.csv")
rownames(BOT2) <- BOT2$Gene.Symbol
BOT2 <- BOT2[,14:16]
BOT2 <- subset(BOT2, rownames(BOT2) %in% DEG_PPI)




#Find the gene names that all datasets have in common
DEG_com <- Reduce(intersect, list(rownames(BOT), rownames(BOT2), rownames(CON),
                                  rownames(DIJ), rownames(DUM), rownames(FFR), 
                                  rownames(LEW),rownames(MID1),rownames(MID2),
                                  rownames(MID3),rownames(MID4),rownames(MOR.FC),
                                  rownames(MOR.SN)))

#Subset each dataset with these common names so they are all the same size
BOT <- subset(BOT, rownames(BOT) %in% DEG_com)
BOT2 <- subset(BOT2, rownames(BOT2) %in% DEG_com)
CON <- subset(CON, rownames(CON) %in% DEG_com)
DIJ <- subset(DIJ, rownames(DIJ) %in% DEG_com)
DUM <- subset(DUM, rownames(DUM) %in% DEG_com)
FFR <- subset(FFR, rownames(FFR) %in% DEG_com)
LEW <- subset(LEW, rownames(LEW) %in% DEG_com)
MID1 <- subset(MID1, rownames(MID1) %in% DEG_com)
MID2 <- subset(MID2, rownames(MID2) %in% DEG_com)
MID3 <- subset(MID3, rownames(MID3) %in% DEG_com)
MID4 <- subset(MID4, rownames(MID4) %in% DEG_com)
MOR.FC <- subset(MOR.FC, rownames(MOR.FC) %in% DEG_com)
MOR.SN <- subset(MOR.SN, rownames(MOR.SN) %in% DEG_com)

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")

write.csv(CON, "CON_DEGfilt.csv")
write.csv(DIJ, "DIJ_DEGfilt.csv")
write.csv(DUM, "DUM_DEGfilt.csv")
write.csv(FFR, "FFR_DEGfilt.csv")
write.csv(LEW, "LEW_DEGfilt.csv")
write.csv(MID1, "MID1_DEGfilt.csv")
write.csv(MID2, "MID2_DEGfilt.csv")
write.csv(MID3, "MID3_DEGfilt.csv")
write.csv(MID4, "MID4_DEGfilt.csv")
write.csv(MOR.FC, "MOR.FC_DEGfilt.csv")
write.csv(MOR.SN, "MOR.SN_DEGfilt.csv")




#### Cor.test Method - Individual files created and analysed on BCBIO ####

library(tictoc)
library(gdata)

##For loop for generating regression values and p values
#### BOT ####
CorExprMat <- t(BOT[1:100,])
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(reg.mat, "BOT_comco.csv")


#### BOT2 ####
CorExprMat <- t(BOT2)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(reg.mat, "BOT2_comco.csv")

#### CON ####
CorExprMat <- t(CON)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "CON_comco.csv")

#### DIJ ####
CorExprMat <- t(DIJ)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "DIJ_comco.csv")

#### DUM ####
CorExprMat <- t(DUM)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "DUM_comco.csv")

#### FFR ####
CorExprMat <- t(FFR)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "FFR_comco.csv")

#### LEW ####
CorExprMat <- t(LEW)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "LEW_comco.csv")

#### MID1 ####
CorExprMat <- t(MID1)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "MID1_comco.csv")

#### MID2 ####
CorExprMat <- t(MID2)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "MID2_comco.csv")

#### MID3 ####
CorExprMat <- t(MID3)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "MID3_comco.csv")

#### MID4 ####
CorExprMat <- t(MID4)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "MID4_comco.csv")


#### MOR.FC ####
CorExprMat <- t(MOR.FC)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "MOR.FC_comco.csv")

#### MOR.SN ####
CorExprMat <- t(MOR.SN)
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")
write.csv(results, "MOR.SN_comco.csv")













#### filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/PD_PPI_BloodOnly/CommonCoexpression/")

#Split names to two columns
CON <- read.csv("CONcorresult.csv")
CON$Gene1 <- as.character(lapply(strsplit(as.character(CON$X), "\\:"), "[", 2))
CON$Gene2 <- as.character(lapply(strsplit(as.character(CON$X), "\\:"), "[", 1))
CON <- CON[,c(5,6,1,2,3,4)]

DIJ <- read.csv("DIJcorresult.csv")
DIJ$Gene1 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 2))
DIJ$Gene2 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 1))
DIJ <- DIJ[,c(5,6,1,2,3,4)]

DUM <- read.csv("DUMcorresult.csv")
DUM$Gene1 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 2))
DUM$Gene2 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 1))
DUM <- DUM[,c(5,6,1,2,3,4)]

FFR <- read.csv("FFRcorresult.csv")
FFR$Gene1 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 2))
FFR$Gene2 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 1))
FFR <- FFR[,c(5,6,1,2,3,4)]
FFR$X <- paste(FFR$Gene1,":",FFR$Gene2, sep = "")

LEW <- read.csv("LEWcorresult.csv")
LEW$Gene1 <- as.character(lapply(strsplit(as.character(LEW$X), "\\:"), "[", 2))
LEW$Gene2 <- as.character(lapply(strsplit(as.character(LEW$X), "\\:"), "[", 1))
LEW <- LEW[,c(5,6,1,2,3,4)]

MID1 <- read.csv("MID1corresult.csv")
MID1$Gene1 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 2))
MID1$Gene2 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 1))
MID1 <- MID1[,c(5,6,1,2,3,4)]
MID1$X <- paste(MID1$Gene1,":",MID1$Gene2, sep = "")

MID2 <- read.csv("MID2corresult.csv")
MID2$Gene1 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 2))
MID2$Gene2 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 1))
MID2 <- MID2[,c(5,6,1,2,3,4)]

MID3 <- read.csv("MID3corresult.csv")
MID3$Gene1 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 2))
MID3$Gene2 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 1))
MID3 <- MID3[,c(5,6,1,2,3,4)]

MID4 <- read.csv("MID4corresult.csv")
MID4$Gene1 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 2))
MID4$Gene2 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 1))
MID4 <- MID4[,c(5,6,1,2,3,4)]

MOR.FC<- read.csv("MOR.FCcorresult.csv")
MOR.FC$Gene1 <- as.character(lapply(strsplit(as.character(MOR.FC$X), "\\:"), "[", 2))
MOR.FC$Gene2 <- as.character(lapply(strsplit(as.character(MOR.FC$X), "\\:"), "[", 1))
MOR.FC <- MOR.FC[,c(5,6,1,2,3,4)]

MOR.SN<- read.csv("MOR.SNcorresult.csv")
MOR.SN$Gene1 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 2))
MOR.SN$Gene2 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 1))
MOR.SN <- MOR.SN[,c(5,6,1,2,3,4)]


### Filter by r value
CON_cor.5 <- CON[CON$reg.mat > 0.5 | CON$reg.mat < -0.5,]
DIJ_cor.5 <- DIJ[DIJ$reg.mat > 0.5 | DIJ$reg.mat < -0.5,]
DUM_cor.5 <- DUM[DUM$reg.mat > 0.5 | DUM$reg.mat < -0.5,]
FFR_cor.5 <- FFR[FFR$reg.mat > 0.5 | FFR$reg.mat < -0.5,]
LEW_cor.5 <- LEW[LEW$reg.mat > 0.5 | LEW$reg.mat < -0.5,]
MID1_cor.5 <- MID1[MID1$reg.mat > 0.5 | MID1$reg.mat < -0.5,]
MID2_cor.5 <- MID2[MID2$reg.mat > 0.5 | MID2$reg.mat < -0.5,]
MID3_cor.5 <- MID3[MID3$reg.mat > 0.5 | MID3$reg.mat < -0.5,]
MID4_cor.5 <- MID4[MID4$reg.mat > 0.5 | MID4$reg.mat < -0.5,]
MOR.FC_cor.5 <- MOR.FC[MOR.FC$reg.mat > 0.5 | MOR.FC$reg.mat < -0.5,]
MOR.SN_cor.5 <- MOR.SN[MOR.SN$reg.mat > 0.5 | MOR.SN$reg.mat < -0.5,]





#Merge into one table. First data is read in, but because the genes have been separated into two columns
#we want them back in one column so that they can be cross referenced correctly. It just so happens that to match
#the cytoscape format, we have to put the gene in the second column first and first column second.


FFR$combo <- paste(FFR$,":",C9_cor.5$Gene1, sep = "")

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
