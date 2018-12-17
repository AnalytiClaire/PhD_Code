#Extract PPI network genes from each dataset
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

##For loop for generating regression values and p values
CorExprMat <- t(VCP[1:100,])
RPTtest <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
p.value <- matrix(0, ncol(test), ncol(test))
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(test)
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    p.value[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$p.value
  }}

rownames(p.value) <- colnames(p.value) <- colnames(test)

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

write.csv(reg.mat, "VCP_coexpression.csv")