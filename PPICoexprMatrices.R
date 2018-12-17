#### Create coexpression matrices from PPI network ####

library(gdata)

C9 <- read.csv("C9uniquegene_samples.csv", row.names = 1)
C9 <- C9[,4:11]
VCP <- read.csv("VCPuniquegene_samples.csv", row.names = 1)
VCP <- VCP[,4:10]
FTLD <- read.csv("FTLDuniquegene_samples.csv", row.names = 1)
FTLD <- FTLD[,9:24]
sALS <- read.csv("sALSuniquegene_samples.csv", row.names = 1)
sALS <- sALS[,4:10]
PET <- read.csv("PET_results_keepfiltering.csv")
rownames(PET) <- PET$hgnc_symbol
PET <- PET[,19:35]
RAV <- read.csv("RAV_results_keepfiltering.csv")
rownames(RAV) <- RAV$hgnc_symbol
RAV <- RAV[,18:30]

#Find the gene names that all datasets have in common
comgene <- read.table("PPI_comgene.txt")

C9.72 <- subset(C9, rownames(C9) %in% comgene)
VCP.72 <- subset(VCP, rownames(VCP) %in% comgene)
FTLD.72 <- subset(FTLD, rownames(FTLD) %in% comgene)
sALS.72 <- subset(sALS, rownames(sALS) %in% comgene)
PET.72 <- subset(PET, rownames(PET) %in% comgene)
RAV.72 <- subset(RAV, rownames(RAV) %in% comgene)
  
###### C9 #########
CorExprMat <- t(C9.72)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}
rownames(reg) <- colnames(reg) <- colnames(test)
ptri <- reg
ptri[lower.tri(reg, diag = TRUE)] <- NA
p.vec <- unmatrix(ptri)
p.vec <- na.omit(p.vec)
p.mat.C9 <- as.data.frame(p.vec)
colnames(p.mat.C9) <- "C9"
p.mat.C9$names <- rownames(p.mat.C9)
p.mat.C9 <- p.mat.C9[order(p.mat.C9$names),]
      

###### sALS #########
CorExprMat <- t(sALS.72)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}
rownames(reg) <- colnames(reg) <- colnames(test)
ptri <- reg
ptri[lower.tri(reg, diag = TRUE)] <- NA
p.vec <- unmatrix(ptri)
p.vec <- na.omit(p.vec)
p.mat.sALS <- as.data.frame(p.vec)
colnames(p.mat.sALS) <- "sALS"
p.mat.sALS$names <- rownames(p.mat.sALS)
p.mat.sALS <- p.mat.sALS[order(p.mat.sALS$names),]

###### FTLD #########
CorExprMat <- t(FTLD.72)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}
rownames(reg) <- colnames(reg) <- colnames(test)
ptri <- reg
ptri[lower.tri(reg, diag = TRUE)] <- NA
p.vec <- unmatrix(ptri)
p.vec <- na.omit(p.vec)
p.mat.FTLD <- as.data.frame(p.vec)
colnames(p.mat.FTLD) <- "FTLD"
p.mat.FTLD$names <- rownames(p.mat.FTLD)
p.mat.FTLD <- p.mat.FTLD[order(p.mat.FTLD$names),]

###### VCP #########
CorExprMat <- t(VCP.72)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}
rownames(reg) <- colnames(reg) <- colnames(test)
ptri <- reg
ptri[lower.tri(reg, diag = TRUE)] <- NA
p.vec <- unmatrix(ptri)
p.vec <- na.omit(p.vec)
p.mat.VCP <- as.data.frame(p.vec)
colnames(p.mat.VCP) <- "VCP"
p.mat.VCP$names <- rownames(p.mat.VCP)
p.mat.VCP <- p.mat.VCP[order(p.mat.VCP$names),]

###### PET #########
CorExprMat <- t(PET.72)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}
rownames(reg) <- colnames(reg) <- colnames(test)
ptri <- reg
ptri[lower.tri(reg, diag = TRUE)] <- NA
p.vec <- unmatrix(ptri)
p.vec <- na.omit(p.vec)
p.mat.PET <- as.data.frame(p.vec)
colnames(p.mat.PET) <- "PET"
p.mat.PET$names <- rownames(p.mat.PET)
p.mat.PET <- p.mat.PET[order(p.mat.PET$names),]

###### RAV #########
CorExprMat <- t(RAV.72)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}
rownames(reg) <- colnames(reg) <- colnames(test)
ptri <- reg
ptri[lower.tri(reg, diag = TRUE)] <- NA
p.vec <- unmatrix(ptri)
p.vec <- na.omit(p.vec)
p.mat.RAV <- as.data.frame(p.vec)
colnames(p.mat.RAV) <- "RAV"
p.mat.RAV$names <- rownames(p.mat.RAV)
p.mat.RAV <- p.mat.RAV[order(p.mat.RAV$names),]


write.csv(p.mat.C9, "pmatC9.csv")
write.csv(p.mat.sALS, "pmatsALS.csv")
write.csv(p.mat.FTLD, "pmatFTLD.csv")
write.csv(p.mat.VCP, "pmatVCP.csv")
write.csv(p.mat.PET, "pmatPET.csv")
write.csv(p.mat.RAV, "pmatRAV.csv")
