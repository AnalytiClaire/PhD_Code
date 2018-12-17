#Selecting DEGS from expression matrix

#Load list of interesting genes
#setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/")
setwd(dir = "/Users/clairegreen/Desktop/")
Genelist <- read.csv("overlap_ens2hgnc_4RNAseq.csv", header = TRUE)

#load dataset
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
exprs <- read.csv("C9rankeduniqueresult.csv")

#Make gene symbol row names
#rownames(exprs) <- exprs$Ensembl
exprspat <- exprs[,49:51]
exprspat[,(length(exprspat)+1)] <- exprs$Ensembl

#Make gene symbol a column
# exprspat <- cbind(exprspat, exprs$Gene.Symbol)
# colnames(exprspat)[length(exprspat)] <- "Gene.Symbol"

#Merge by interesting gene names with expression to form matrix
patgene <- merge(Genelist, exprspat, by.x = "ensembl_gene_id", by.y = "V4")
#patgene <- patgene[!duplicated(patgene[,11]),]
# rownames(patgene) <- patgene$V1
# patgene[,1] <- NULL

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/Co-expression/")
write.csv(patgene, file = "VCP_DEG_CON_ens.csv")

#### Cor.test Method ####

library(tictoc)
library(gdata)

##load dataset
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/Co-expression/")
Exprs_val <- read.csv("C9r_DEG_CON_Exprs.csv")
rownames(Exprs_val) <- Exprs_val[,2]
Exprs_val[,1:2] <- NULL

CorExprMat <- t(Exprs_val)

reg <- matrix(0, ncol(CorExprMat), ncol(CorExprMat))

tic()
for (i in 1:ncol(CorExprMat)){
  for (j in 1:ncol(CorExprMat)){
    reg[i,j] <- cor.test(CorExprMat[,i], CorExprMat[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(CorExprMat)
toc()


#Extract R values
corRadj <- reg
corRadj[lower.tri(corRadj, diag = TRUE)] <- NA

#Turn into vector
corRadj <- as.matrix(corRadj)
corRvec <- unmatrix(corRadj)
#Remove NA values
corRvec <- na.omit(corRvec)
corRvec <- as.data.frame(corRvec)

write.csv(corRvec, file = "C9r_CON_CorResults.csv")


#Generate matrix with all Rho values

C9mR <- read.csv("C9m_CorResults.csv")
CHmR <- read.csv("CHMP2B_CorResults.csv")
GRNR <- read.csv("GRN_CorResults.csv")
VCPR <- read.csv("VCP_CorResults.csv")
C9rR <- read.csv("C9r_CorResults.csv")

C9mCON <- read.csv("C9m_CON_CorResults.csv")
C9rCON <- read.csv("C9r_CON_CorResults.csv")
CHCON <- read.csv("CH_CON_CorResults.csv")
GRNCON <- read.csv("FTLD_CON_CorResults.csv")
VCPCON <- read.csv("VCP_CON_CorResults.csv")



RhoValues <- merge(C9mR, CHmR, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, GRNR, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, VCPR, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, C9rR, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, C9mCON, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, C9rCON, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, CHCON, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, GRNCON, by.x = "X", by.y = "X")
RhoValues <- merge(RhoValues, VCPCON, by.x = "X", by.y = "X")

colnames(RhoValues) <- c("GenePair","C9m", "CHMP2B", "GRN", "VCP", "C9r", "C9mCON", "C9rCON",
                         "CHCON", "GRNCON", "VCPCON")

rownames(RhoValues) <- RhoValues$GenePair
RhoValues[,1] <- NULL

Rho <- matrix(0, ncol(RhoValues), ncol(RhoValues))

tic()
for (i in 1:ncol(RhoValues)){
  for (j in 1:ncol(RhoValues)){
    Rho[i,j] <- cor.test(RhoValues[,i], RhoValues[,j], method = "kendall")$p.value
  }}

rownames(Rho) <- colnames(Rho) <- colnames(RhoValues)
toc()


#Extract R values
corRadj <- Rho
corRadj[lower.tri(corRadj, diag = TRUE)] <- NA

#Turn into vector
corRadj <- as.matrix(corRadj)
corRvec <- unmatrix(corRadj)
#Remove NA values
corRvec <- na.omit(corRvec)
corRvec <- as.data.frame(corRvec)

write.csv(corRvec, file = "Alldatasets_Pval.csv")





## Rank values by Spearman's Rho
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/Co-expression/")
# 
# C9mR <- read.csv("C9m_CorResults.csv")
# rankC9m <- C9mR[order(C9mR[,1], decreasing = FALSE),]
# #rankC9m$rank <- seq.int(nrow(rankC9m))
# 
# CHmR <- read.csv("CHMP2B_CorResults.csv")
# rankCH <- CHmR[order(CHmR[,1], decreasing = FALSE),]
# #rankCH$rank <- seq.int(nrow(rankCH))
# 
# GRNR <- read.csv("GRN_CorResults.csv")
# rankGRN <- GRNR[order(GRNR[,1], decreasing = FALSE),]
# #rankGRN$rank <- seq.int(nrow(rankGRN))
# 
# VCPR <- read.csv("VCP_CorResults.csv")
# rankVCP <- VCPR[order(VCPR[,1], decreasing = FALSE),]
# #rankVCP$rank <- seq.int(nrow(rankVCP))
# 
# C9rR <- read.csv("C9r_CorResults.csv")
# rankC9r <- C9rR[order(C9rR[,1], decreasing = FALSE),]
# #rankC9r$rank <- seq.int(nrow(rankC9r))
# 
# 
# RankC9m_CHm <- cor.test(rankC9m$corRvec, rankCH$corRvec, method = "kendall")
# RankC9m_GRN <- cor.test(rankC9m$corRvec, rankGRN$corRvec, method = "kendall")
# RankC9m_VCP <- cor.test(rankC9m$corRvec, rankVCP$corRvec, method = "kendall")
# RankC9m_C9r <- cor.test(rankC9m$corRvec, rankC9r$corRvec, method = "kendall")
# RankCHm_GRN <- cor.test(rankCH$corRvec, rankGRN$corRvec, method = "kendall")
# RankCHm_VCP <- cor.test(rankCH$corRvec, rankVCP$corRvec, method = "kendall")
# RankCHm_C9r <- cor.test(rankCH$corRvec, rankC9r$corRvec, method = "kendall")
# RankGRN_VCP <- cor.test(rankGRN$corRvec, rankVCP$corRvec, method = "kendall")
# RankGRN_C9r <- cor.test(rankGRN$corRvec, rankC9r$corRvec, method = "kendall")
# RankVCP_C9r <- cor.test(rankVCP$corRvec, rankC9r$corRvec, method = "kendall")
# 
# Taus <- cbind(RankC9m_CHm$estimate, RankC9m_GRN$estimate, RankC9m_VCP$estimate, 
#               RankC9m_C9r$estimate, RankCHm_GRN$estimate, RankCHm_VCP$estimate, 
#               RankCHm_C9r$estimate, RankGRN_VCP$estimate, RankGRN_C9r$estimate, 
#               RankVCP_C9r$estimate)
# colnames(Taus) <- c("RankC9m_CHm", "RankC9m_GRN", "RankC9m_VCP", "RankC9m_C9r",
#                     "RankCHm_GRN", "RankCHm_VCP", "RankCHm_C9r", "RankGRN_VCP", 
#                     "RankGRN_C9r", "RankVCP_C9r")
# 
# Pval <- cbind(RankC9m_CHm$p.value, RankC9m_GRN$p.value, RankC9m_VCP$p.value, 
#               RankC9m_C9r$p.value, RankCHm_GRN$p.value, RankCHm_VCP$p.value, 
#               RankCHm_C9r$p.value, RankGRN_VCP$p.value, RankGRN_C9r$p.value, 
#               RankVCP_C9r$p.value)
# 
# colnames(Pval) <- c("RankC9m_CHm", "RankC9m_GRN", "RankC9m_VCP", "RankC9m_C9r",
#                     "RankCHm_GRN", "RankCHm_VCP", "RankCHm_C9r", "RankGRN_VCP", 
#                     "RankGRN_C9r", "RankVCP_C9r")
# 
# result <- rbind(Taus, Pval)
# row.names(result) <- c("Tau", "pValue")
# 
# result <- t(result)
# 
# write.csv(result, "DatasetCorrelationResults_Kendall.csv")
