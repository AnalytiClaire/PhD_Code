#### Parkinson's DEG PPI Correlation ####

#Read in network nodes

#### Familial Blood ####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PPIGenes.txt")

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
#Extract PPI network genes from each dataset#####
LEW <- read.csv("LEWfilteredresult.csv")
rownames(LEW) <- LEW$Gene.Symbol
LEW <- LEW[,24:29]
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




#Find the gene names that all datasets have in common
DEG_com <- Reduce(intersect, list(rownames(DIJ), rownames(FFR), 
                                  rownames(LEW),rownames(MID1),rownames(MID2),
                                  rownames(MID3),rownames(MID4),rownames(MOR.FC),
                                  rownames(MOR.SN), rownames(DUM)))

#Subset each dataset with these common names so they are all the same size
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

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/")

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

#Run Correlation analysis on Sharc

####
###TEST####
# uniqueresult <- read.csv("DIJ_DEGfilt.csv", row.names = 1)

# ##For loop for generating regression values and p values
# CorExprMat <- t(uniqueresult)
# 
# reg <- matrix(0, ncol(CorExprMat), ncol(CorExprMat))
# p.value <- matrix(0, ncol(CorExprMat), ncol(CorExprMat))
# 
# for (i in 1:ncol(CorExprMat)){
#   for (j in 1:ncol(CorExprMat)){
#     reg[i,j] <- cor.test(CorExprMat[,i], CorExprMat[,j], method = "spearman")$estimate
#   }}
# 
# rownames(reg) <- colnames(reg) <- colnames(CorExprMat)
# 
# for (i in 1:ncol(CorExprMat)){
#   for (j in 1:ncol(CorExprMat)){
#     p.value[i,j] <- cor.test(CorExprMat[,i], CorExprMat[,j], method = "pearson")$p.value
#   }}
# 
# rownames(p.value) <- colnames(p.value) <- colnames(CorExprMat)
# 
# ##Only take upper triangle without diagonal (all comparisons are currently doubled)
# ptri <- p.value
# ptri[lower.tri(ptri, diag = TRUE)] <- NA
# 
# #Turn into vector
# library(gdata)
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
# write.csv(results, "testDIJcorresult.csv")
####

#### filter correlations
setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/")

DIJ <- read.csv("testDIJcorresult.csv")
DIJ$Gene1 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 2))
DIJ$Gene2 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 1))
DIJ <- DIJ[,c(5,6,1,2,3,4)]
# DIJ$X <- paste(DIJ$Gene1,":",DIJ$Gene2, sep = "")

DUM <- read.csv("DUMcorresult.csv")
DUM$Gene1 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 2))
DUM$Gene2 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 1))
DUM <- DUM[,c(5,6,1,2,3,4)]
DUM$X <- paste(DUM$Gene1,":",DUM$Gene2, sep = "")

FFR <- read.csv("FFRcorresult.csv")
FFR$Gene1 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 2))
FFR$Gene2 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 1))
FFR <- FFR[,c(5,6,1,2,3,4)]
# FFR$X <- paste(FFR$Gene1,":",FFR$Gene2, sep = "")

LEW <- read.csv("LEWcorresult.csv")
LEW$Gene1 <- as.character(lapply(strsplit(as.character(LEW$X), "\\:"), "[", 2))
LEW$Gene2 <- as.character(lapply(strsplit(as.character(LEW$X), "\\:"), "[", 1))
LEW <- LEW[,c(5,6,1,2,3,4)]
LEW$X <- paste(LEW$Gene1,":",LEW$Gene2, sep = "")

MID1 <- read.csv("MID1corresult.csv")
MID1$Gene1 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 2))
MID1$Gene2 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 1))
MID1 <- MID1[,c(5,6,1,2,3,4)]
# MID1$X <- paste(MID1$Gene1,":",MID1$Gene2, sep = "")
# 
# MID2 <- read.csv("MID2corresult.csv")
# MID2$Gene1 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 2))
# MID2$Gene2 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 1))
# MID2 <- MID2[,c(5,6,1,2,3,4)]
# MID2$X <- paste(MID2$Gene2,":",MID2$Gene1, sep = "")
# 
# MID3 <- read.csv("MID3corresult.csv")
# MID3$Gene1 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 2))
# MID3$Gene2 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 1))
# MID3 <- MID3[,c(5,6,1,2,3,4)]
# # MID3$X <- paste(MID3$Gene1,":",MID3$Gene2, sep = "")
# 
# 
# MID4 <- read.csv("MID4corresult.csv")
# MID4$Gene1 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 2))
# MID4$Gene2 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 1))
# MID4 <- MID4[,c(5,6,1,2,3,4)]
# # MID4$X <- paste(MID4$Gene1,":",MID4$Gene2, sep = "")
# 
# MOR.FC<- read.csv("MOR.FCcorresult.csv")
# MOR.FC$Gene1 <- as.character(lapply(strsplit(as.character(MOR.FC$X), "\\:"), "[", 2))
# MOR.FC$Gene2 <- as.character(lapply(strsplit(as.character(MOR.FC$X), "\\:"), "[", 1))
# MOR.FC <- MOR.FC[,c(5,6,1,2,3,4)]
# # MOR.FC$X <- paste(MOR.FC$Gene1,":",MOR.FC$Gene2, sep = "")
# 
# MOR.SN<- read.csv("MOR.SNcorresult.csv")
# MOR.SN$Gene1 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 2))
# MOR.SN$Gene2 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 1))
# MOR.SN <- MOR.SN[,c(5,6,1,2,3,4)]
# # MOR.SN$X <- paste(MOR.SN$Gene1,":",MOR.SN$Gene2, sep = "")



# 
# write.csv(DIJ, "DIJ_coexpr.csv")
# write.csv(DUM, "DUM_coexpr.csv")
# write.csv(FFR, "FFR_coexpr.csv")
# write.csv(LEW, "LEW_coexpr.csv")
# write.csv(MID1, "MID1_coexpr.csv")
# write.csv(MID2, "MID2_coexpr.csv")
# write.csv(MID3, "MID3_coexpr.csv")
# write.csv(MID4, "MID4_coexpr.csv")
# write.csv(MOR.FC, "MOR.FC_coexpr.csv")
# write.csv(MOR.SN, "MOR.SN_coexpr.csv")


setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/CorResultSepGeneName/")

DIJ <- read.csv("DIJ_coexpr.csv", row.names = 1)
DUM <- read.csv("DUM_coexpr.csv", row.names = 1)
FFR <- read.csv("FFR_coexpr.csv", row.names = 1)
LEW <- read.csv("LEW_coexpr.csv", row.names = 1)
MID1 <- read.csv("MID1_coexpr.csv", row.names = 1)
MID2 <- read.csv("MID2_coexpr.csv", row.names = 1)
MID3 <- read.csv("MID3_coexpr.csv", row.names = 1)
MID4 <- read.csv("MID4_coexpr.csv", row.names = 1)
MOR.FC <- read.csv("MOR.FC_coexpr.csv", row.names = 1)
MOR.SN <- read.csv("MOR.SN_coexpr.csv", row.names = 1)





hist(DIJ$reg.mat, main = "Coexpression Distribution (DIJ)", xlab = "Coexpression Value")
hist(DUM$reg.mat, main = "Coexpression Distribution (DUM)", xlab = "Coexpression Value")
hist(FFR$reg.mat, main = "Coexpression Distribution (FFR)", xlab = "Coexpression Value")
hist(LEW$reg.mat, main = "Coexpression Distribution (LEW)", xlab = "Coexpression Value")
hist(MID1$reg.mat, main = "Coexpression Distribution (MID1)", xlab = "Coexpression Value")
hist(MID2$reg.mat, main = "Coexpression Distribution (MID2)", xlab = "Coexpression Value") #skew
hist(MID3$reg.mat, main = "Coexpression Distribution (MID3)", xlab = "Coexpression Value")
hist(MID4$reg.mat, main = "Coexpression Distribution (MID4)", xlab = "Coexpression Value")
hist(MOR.FC$reg.mat, main = "Coexpression Distribution (MOR.FC)", xlab = "Coexpression Value")
hist(MOR.SN$reg.mat, main = "Coexpression Distribution (MOR.SN)", xlab = "Coexpression Value")




thresh <- 0.5
### Filter by r value
DIJ_cor.5 <- DIJ[DIJ$reg.mat > thresh | DIJ$reg.mat < -thresh,]
DUM_cor.5 <- DUM[DUM$reg.mat > thresh | DUM$reg.mat < -thresh,]
FFR_cor.5 <- FFR[FFR$reg.mat > thresh | FFR$reg.mat < -thresh,]
# LEW_cor.5 <- LEW[LEW$reg.mat > thresh | LEW$reg.mat < -thresh,]
MID1_cor.5 <- MID1[MID1$reg.mat > thresh | MID1$reg.mat < -thresh,]
MID2_cor.5 <- MID2[MID2$reg.mat > thresh | MID2$reg.mat < -thresh,]
MID3_cor.5 <- MID3[MID3$reg.mat > thresh | MID3$reg.mat < -thresh,]
MID4_cor.5 <- MID4[MID4$reg.mat > thresh | MID4$reg.mat < -thresh,]
# MOR.FC_cor.5 <- MOR.FC[MOR.FC$reg.mat > thresh | MOR.FC$reg.mat < -thresh,]
MOR.SN_cor.5 <- MOR.SN[MOR.SN$reg.mat > thresh | MOR.SN$reg.mat < -thresh,]



DIJ_cor.5$X2 <- paste(DIJ_cor.5$Gene2,":",DIJ_cor.5$Gene1, sep = "")
DUM_cor.5$X2 <- paste(DUM_cor.5$Gene2,":",DUM_cor.5$Gene1, sep = "")
FFR_cor.5$X2 <- paste(FFR_cor.5$Gene2,":",FFR_cor.5$Gene1, sep = "")
MID1_cor.5$X2 <- paste(MID1_cor.5$Gene1,":",MID1_cor.5$Gene2, sep = "")
MID2_cor.5$X2 <- paste(MID2_cor.5$Gene1,":",MID2_cor.5$Gene2, sep = "")
MID3_cor.5$X2 <- paste(MID3_cor.5$Gene1,":",MID3_cor.5$Gene2, sep = "")
MID4_cor.5$X2 <- paste(MID4_cor.5$Gene2,":",MID4_cor.5$Gene1, sep = "")
MOR.SN_cor.5$X2 <- paste(MOR.SN_cor.5$Gene1,":",MOR.SN_cor.5$Gene2, sep = "")

#Merge first two datasets
corresult1 <- merge(DIJ_cor.5, FFR_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(DIJ_cor.5, FFR_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR")

#Merge with 3rd Dataset
corresult3 <- merge(result, MID1_cor.5, by.x = "Correlation", by.y = "X")
corresult4 <- merge(result, MID1_cor.5, by.x = "Correlation", by.y = "X2")

corresult3 <- corresult3[, c(1:5, 10)]
corresult4 <- corresult4[, c(1:5, 11)]

result <- rbind(corresult3, corresult4)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1")

#Merge with 4th Dataset
corresult5 <- merge(result, MID2_cor.5, by.x = "Correlation", by.y = "X")
corresult6 <- merge(result, MID2_cor.5, by.x = "Correlation", by.y = "X2")

corresult5 <- corresult5[, c(1:6, 11)]
corresult6 <- corresult6[, c(1:6, 12)]

result <- rbind(corresult5, corresult6)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2")

#Merge with 5th Dataset
corresult7 <- merge(result, MID3_cor.5, by.x = "Correlation", by.y = "X")
corresult8 <- merge(result, MID3_cor.5, by.x = "Correlation", by.y = "X2")

corresult7 <- corresult7[, c(1:7, 12)]
corresult8 <- corresult8[, c(1:7, 13)]

result <- rbind(corresult7, corresult8)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3")

#Merge with 6th Dataset
corresult9 <- merge(result, MID4_cor.5, by.x = "Correlation", by.y = "X")
corresult10 <- merge(result, MID4_cor.5, by.x = "Correlation", by.y = "X2")

corresult9 <- corresult9[, c(1:8, 13)]
corresult10 <- corresult10[, c(1:8, 14)]

result <- rbind(corresult9, corresult10)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3", "MID4")

#Merge with 7th Dataset
corresult11 <- merge(result, MOR.SN_cor.5, by.x = "Correlation", by.y = "X")
corresult12 <- merge(result, MOR.SN_cor.5, by.x = "Correlation", by.y = "X2")

corresult11 <- corresult11[, c(1:9, 14)]
corresult12 <- corresult12[, c(1:9, 15)]

result <- rbind(corresult11, corresult12)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3", "MID4", "MOR.SN")

#Merge with 8th Dataset
corresult13 <- merge(result, DUM_cor.5, by.x = "Correlation", by.y = "X")
corresult14 <- merge(result, DUM_cor.5, by.x = "Correlation", by.y = "X2")

corresult13 <- corresult13[, c(1:10, 15)]
corresult14 <- corresult14[, c(1:10, 16)]

result <- rbind(corresult13, corresult14)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3", "MID4", "MOR.SN")






# 
# 
# write.csv(DIJ_cor.5, "DIJ_coexpr_point3.csv")
# write.csv(FFR_cor.5, "FFR_coexpr_point3.csv")
# write.csv(MID1_cor.5, "MID1_coexpr_point3.csv")
# write.csv(MID2_cor.5, "MID2_coexpr_point3.csv")
# write.csv(MID3_cor.5, "MID3_coexpr_point3.csv")
# write.csv(MID4_cor.5, "MID4_coexpr_point3.csv")
# write.csv(MOR.SN_cor.5, "MOR.SN_coexpr_point3.csv")
# 
# 
# 
# ### Find matches 
# 
# Commonedge <- Reduce(intersect, list(DIJ_cor.5$X, FFR_cor.5$X, MID1_cor.5$X, 
#                                      MID2_cor.5$X,MID3_cor.5$X, MID4_cor.5$X,
#                                      MOR.SN_cor.5$X))
# 
# 
# #Subset each dataset with these common names so they are all the same size
# DIJ_CE <- subset(DIJ_cor.5, DIJ_cor.5$X %in% Commonedge)
# DIJ_CE <- DIJ_CE[order(DIJ_CE$X),]
# # DUM_CE <- subset(DUM_cor.5, DUM_cor.5$X %in% Commonedge)
# # DUM_CE <- DUM_CE[order(DUM_CE$X),]
# FFR_CE <- subset(FFR_cor.5, FFR_cor.5$X %in% Commonedge)
# FFR_CE <- FFR_CE[order(FFR_CE$X),]
# # LEW_CE <- subset(LEW_cor.5, LEW_cor.5$X %in% Commonedge)
# # LEW_CE <- LEW_CE[order(LEW_CE$X),]
# MID1_CE <- subset(MID1_cor.5, MID1_cor.5$X %in% Commonedge)
# MID1_CE <- MID1_CE[order(MID1_CE$X),]
# MID2_CE <- subset(MID2_cor.5, MID2_cor.5$X %in% Commonedge)
# MID2_CE <- MID2_CE[order(MID2_CE$X),]
# MID3_CE <- subset(MID3_cor.5, MID3_cor.5$X %in% Commonedge)
# MID3_CE <- MID3_CE[order(MID3_CE$X),]
# MID4_CE <- subset(MID4_cor.5, MID4_cor.5$X %in% Commonedge)
# MID4_CE <- MID4_CE[order(MID4_CE$X),]
# # MOR.FC_CE <- subset(MOR.FC_cor.5, MOR.FC_cor.5$X %in% Commonedge)
# # MOR.FC_CE <- MOR.FC_CE[order(MOR.FC_CE$X),]
# MOR.SN_CE <- subset(MOR.SN_cor.5, MOR.SN_cor.5$X %in% Commonedge)
# MOR.SN_CE <- MOR.SN_CE[order(MOR.SN_CE$X),]
# 
# 
# 
# CommonGroup <- data.frame(row.names = DIJ_CE$X,
#                           DIJ = DIJ_CE$reg.mat,
#                           FFR = FFR_CE$reg.mat,
#                           MID1 = MID1_CE$reg.mat,
#                           MID2 = MID2_CE$reg.mat, 
#                           MID3 = MID3_CE$reg.mat, 
#                           MID4 = MID4_CE$reg.mat, 
#                           MOR.SN = MOR.SN_CE$reg.mat)

CommonGroup <- result
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:10]

CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "PD_0point5_wDUM.csv", quote = F, row.names = F)


x <- readLines("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PDNuggetGenes.txt")
y <- readLines("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/0point3nuggetgenes.txt")


length(intersect(x,y))


nodetable <- read.csv("PD_samedir_1node.csv")
nuggenes <- as.character(nodetable$name)
PDgenes <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/ParkinsonsDiseaseMalacards.txt")
celltype <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv")
DEG <- readLines("ALS_sfblood_ALLgenes.txt")
nalls <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/NallsPDGWAS.txt")
targetvalPD <- c("SIRT2", "DNM1L", "STMN1", "DNAJA3", "TBK1", "RTCA", "ANXA7", "DNAJC12", "RTN4", "ADSL", "MDH1","ATP6V1G2", 
                 "YWHAZ", "ETS1")

nodetable$PDMalacards <- nodetable$shared.name %in% PDgenes
nodetable$DEG <- nodetable$name %in% DEG
nodetable$targetvalPD <- nodetable$shared.name %in% targetvalPD
nodetable$targetvalPD <- nodetable$shared.name %in% nalls

nodetable_celltype <- merge(celltype, nodetable, by.x = "Gene.symbol",  by.y = "shared.name", all = T)
nodetable_celltype <- subset(nodetable_celltype, !(nodetable_celltype$SUID == "NA"))


nodetable_celltype$targetvalPD <- nodetable_celltype$Gene.symbol %in% targetvalPD
write.csv(nodetable_celltype, "ModifiedPDNodetable.csv", row.names = F)

### Control network ####
#### Familial Blood ####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PDNuggetGenes.txt")

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

DIJ <- read.csv("DIJfilteredresult.csv")
rownames(DIJ) <- DIJ$Gene.Symbol
DIJ <- DIJ[,21:28]
DIJ <- subset(DIJ, rownames(DIJ) %in% DEG_PPI)

DUM <- read.csv("DUM_UniqueGene_DESeq2.csv")
rownames(DUM) <- DUM$hgnc_symbol
DUM <- DUM[,9:52]
DUM <- subset(DUM, rownames(DUM) %in% DEG_PPI)

FFR <- read.csv("FFRfilteredresult.csv")
rownames(FFR) <- FFR$Gene.Symbol
FFR <- FFR[,20:28]
FFR <- subset(FFR, rownames(FFR) %in% DEG_PPI)

MID1 <- read.csv("MID1filteredresult.csv")
rownames(MID1) <- MID1$Gene.Symbol
MID1 <- MID1[,20:27]
MID1 <- subset(MID1, rownames(MID1) %in% DEG_PPI)

MID2 <- read.csv("MID2filteredresult.csv")
rownames(MID2) <- MID2$Gene.Symbol
MID2 <- MID2[,21:38]
MID2 <- subset(MID2, rownames(MID2) %in% DEG_PPI)

MID3 <- read.csv("MID3filteredresult.csv")
rownames(MID3) <- MID3$Gene.Symbol
MID3 <- MID3[,21:35]
MID3 <- subset(MID3, rownames(MID3) %in% DEG_PPI)

MID4 <- read.csv("MID4filteredresult.csv")
rownames(MID4) <- MID4$Gene.Symbol
MID4 <- MID4[,21:40]
MID4 <- subset(MID4, rownames(MID4) %in% DEG_PPI)

MOR.SN <- read.csv("MOR.SNfilteredresult.csv")
rownames(MOR.SN) <- MOR.SN$Gene.Symbol
MOR.SN <- MOR.SN[,21:35]
MOR.SN <- subset(MOR.SN, rownames(MOR.SN) %in% DEG_PPI)

setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression")

write.csv(DIJ, "DIJ_DEGfilt.csv")
write.csv(FFR, "FFR_DEGfilt.csv")
write.csv(DUM, "DUM_DEGfilt.csv")
write.csv(MID1, "MID1_DEGfilt.csv")
write.csv(MID2, "MID2_DEGfilt.csv")
write.csv(MID3, "MID3_DEGfilt.csv")
write.csv(MID4, "MID4_DEGfilt.csv")
write.csv(MOR.SN, "MOR.SN_DEGfilt.csv")



#### Cor.test Method ####

library(tictoc)
library(gdata)

##For loop for generating regression values and p values
CorExprMat <- t(MOR.SN)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
p.value <- matrix(0, ncol(test), ncol(test))

Sys.time()
tic()
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(test)
toc()
Sys.time()
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
write.csv(results, "MOR.SN_con_coexpression.csv")


#### filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")

DIJ <- read.csv("DIJ_con_coexpression.csv")
DIJ$Gene1 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 2))
DIJ$Gene2 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 1))
DIJ <- DIJ[,c(1,5,6,2,3,4)]

FFR <- read.csv("FFR_con_coexpression.csv")
FFR$Gene1 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 2))
FFR$Gene2 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 1))
FFR <- FFR[,c(1,5,6,2,3,4)]

DUM <- read.csv("DUM_con_coexpression.csv")
DUM$Gene1 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 2))
DUM$Gene2 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 1))
DUM <- DUM[,c(1,5,6,2,3,4)]

MID1 <- read.csv("MID1_con_coexpression.csv")
MID1$Gene1 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 2))
MID1$Gene2 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 1))
MID1 <- MID1[,c(1,5,6,2,3,4)]

MID2 <- read.csv("MID2_con_coexpression.csv")
MID2$Gene1 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 2))
MID2$Gene2 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 1))
MID2 <- MID2[,c(1,5,6,2,3,4)]

MID3 <- read.csv("MID3_con_coexpression.csv")
MID3$Gene1 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 2))
MID3$Gene2 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 1))
MID3 <- MID3[,c(1,5,6,2,3,4)]

MID4 <- read.csv("MID4_con_coexpression.csv")
MID4$Gene1 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 2))
MID4$Gene2 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 1))
MID4 <- MID4[,c(1,5,6,2,3,4)]

MOR.SN<- read.csv("MOR.SN_con_coexpression.csv")
MOR.SN$Gene1 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 2))
MOR.SN$Gene2 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 1))
MOR.SN <- MOR.SN[,c(1,5,6,2,3,4)]



thresh <- 0.5
### Filter by r value
DIJ_cor.5 <- DIJ[DIJ$reg.mat > thresh | DIJ$reg.mat < -thresh,]
DUM_cor.5 <- DUM[DUM$reg.mat > thresh | DUM$reg.mat < -thresh,]
FFR_cor.5 <- FFR[FFR$reg.mat > thresh | FFR$reg.mat < -thresh,]
MID1_cor.5 <- MID1[MID1$reg.mat > thresh | MID1$reg.mat < -thresh,]
MID2_cor.5 <- MID2[MID2$reg.mat > thresh | MID2$reg.mat < -thresh,]
MID3_cor.5 <- MID3[MID3$reg.mat > thresh | MID3$reg.mat < -thresh,]
MID4_cor.5 <- MID4[MID4$reg.mat > thresh | MID4$reg.mat < -thresh,]
MOR.SN_cor.5 <- MOR.SN[MOR.SN$reg.mat > thresh | MOR.SN$reg.mat < -thresh,]


DIJ_cor.5$X2 <- paste(DIJ_cor.5$Gene1,":",DIJ_cor.5$Gene2, sep = "")
DUM_cor.5$X2 <- paste(DUM_cor.5$Gene1,":",DUM_cor.5$Gene2, sep = "")
FFR_cor.5$X2 <- paste(FFR_cor.5$Gene1,":",FFR_cor.5$Gene2, sep = "")
MID1_cor.5$X2 <- paste(MID1_cor.5$Gene1,":",MID1_cor.5$Gene2, sep = "")
MID2_cor.5$X2 <- paste(MID2_cor.5$Gene1,":",MID2_cor.5$Gene2, sep = "")
MID3_cor.5$X2 <- paste(MID3_cor.5$Gene1,":",MID3_cor.5$Gene2, sep = "")
MID4_cor.5$X2 <- paste(MID4_cor.5$Gene1,":",MID4_cor.5$Gene2, sep = "")
MOR.SN_cor.5$X2 <- paste(MOR.SN_cor.5$Gene1,":",MOR.SN_cor.5$Gene2, sep = "")

#Merge first two datasets
corresult1 <- merge(DIJ_cor.5, FFR_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(DIJ_cor.5, FFR_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(2:4, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR")

#Merge with 3rd Dataset
corresult3 <- merge(result, MID1_cor.5, by.x = "Correlation", by.y = "X")
corresult4 <- merge(result, MID1_cor.5, by.x = "Correlation", by.y = "X2")

corresult3 <- corresult3[, c(1:5, 10)]
corresult4 <- corresult4[, c(1:5, 11)]

result <- rbind(corresult3, corresult4)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1")

#Merge with 4th Dataset
corresult5 <- merge(result, MID2_cor.5, by.x = "Correlation", by.y = "X")
corresult6 <- merge(result, MID2_cor.5, by.x = "Correlation", by.y = "X2")

corresult5 <- corresult5[, c(1:6, 11)]
corresult6 <- corresult6[, c(1:6, 12)]

result <- rbind(corresult5, corresult6)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2")

#Merge with 5th Dataset
corresult7 <- merge(result, MID3_cor.5, by.x = "Correlation", by.y = "X")
corresult8 <- merge(result, MID3_cor.5, by.x = "Correlation", by.y = "X2")

corresult7 <- corresult7[, c(1:7, 12)]
corresult8 <- corresult8[, c(1:7, 13)]

result <- rbind(corresult7, corresult8)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3")

#Merge with 6th Dataset
corresult9 <- merge(result, MID4_cor.5, by.x = "Correlation", by.y = "X")
corresult10 <- merge(result, MID4_cor.5, by.x = "Correlation", by.y = "X2")

corresult9 <- corresult9[, c(1:8, 13)]
corresult10 <- corresult10[, c(1:8, 14)]

result <- rbind(corresult9, corresult10)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3", "MID4")

#Merge with 7th Dataset
corresult11 <- merge(result, MOR.SN_cor.5, by.x = "Correlation", by.y = "X")
corresult12 <- merge(result, MOR.SN_cor.5, by.x = "Correlation", by.y = "X2")

corresult11 <- corresult11[, c(1:9, 14)]
corresult12 <- corresult12[, c(1:9, 15)]

result <- rbind(corresult11, corresult12)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3", "MID4", "MOR.SN")

#Merge with 8th Dataset
corresult13 <- merge(result, DUM_cor.5, by.x = "Correlation", by.y = "X")
corresult14 <- merge(result, DUM_cor.5, by.x = "Correlation", by.y = "X2")

corresult13 <- corresult13[, c(1:10, 15)]
corresult14 <- corresult14[, c(1:10, 16)]

result <- rbind(corresult13, corresult14)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "DIJ", "FFR", "MID1", "MID2", "MID3", "MID4", "MOR.SN", "DUM")
# 
# 
# 
# ### Find matches 
# 
# Commonedge <- Reduce(intersect, list(DIJ_cor.5$X,FFR_cor.5$X, MID1_cor.5$X, 
#                                      MID2_cor.5$X,MID3_cor.5$X, MID4_cor.5$X,
#                                      MOR.SN_cor.5$X))
# 
# 
# #Subset each dataset with these common names so they are all the same size
# DIJ_CE <- subset(DIJ_cor.5, DIJ_cor.5$X %in% Commonedge)
# DIJ_CE <- DIJ_CE[order(DIJ_CE$X),]
# FFR_CE <- subset(FFR_cor.5, FFR_cor.5$X %in% Commonedge)
# FFR_CE <- FFR_CE[order(FFR_CE$X),]
# MID1_CE <- subset(MID1_cor.5, MID1_cor.5$X %in% Commonedge)
# MID1_CE <- MID1_CE[order(MID1_CE$X),]
# MID2_CE <- subset(MID2_cor.5, MID2_cor.5$X %in% Commonedge)
# MID2_CE <- MID2_CE[order(MID2_CE$X),]
# MID3_CE <- subset(MID3_cor.5, MID3_cor.5$X %in% Commonedge)
# MID3_CE <- MID3_CE[order(MID3_CE$X),]
# MID4_CE <- subset(MID4_cor.5, MID4_cor.5$X %in% Commonedge)
# MID4_CE <- MID4_CE[order(MID4_CE$X),]
# MOR.SN_CE <- subset(MOR.SN_cor.5, MOR.SN_cor.5$X %in% Commonedge)
# MOR.SN_CE <- MOR.SN_CE[order(MOR.SN_CE$X),]



CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:11]

CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "Control_0point5_2.csv", quote = F, row.names = F)


pat <- read.csv("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PD_0point5_wDUM.csv")
con <- read.csv("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/Control_0point5_2.csv")

patcon <- merge(pat, con, by = "Gene")
patcon <- patcon[,c(1,9,20)]
patcon$Gene1 <- as.character(lapply(strsplit(as.character(patcon$Gene), "\\:"), "[", 2))
patcon$Gene2 <- as.character(lapply(strsplit(as.character(patcon$Gene), "\\:"), "[", 1))
write.csv(patcon, "overlap_remove.csv", row.names = F, quote = F)

#### Familial Blood Network ####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PD_NuggetGenes_conoverlapremoved.txt")

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

ATP13A2 <- read.csv("AMA_ATP13A2rankeduniqueresult.csv")
rownames(ATP13A2) <- ATP13A2$Gene.Symbol
ATP13A2 <- ATP13A2[,204:208]
ATP13A2 <- subset(ATP13A2, rownames(ATP13A2) %in% DEG_PPI)

PARKIN <- read.csv("AMA_PARKINrankeduniqueresult.csv")
rownames(PARKIN) <- PARKIN$Gene.Symbol
PARKIN <- PARKIN[,204:216]
PARKIN <- subset(PARKIN, rownames(PARKIN) %in% DEG_PPI)

PINK1 <- read.csv("AMA_PINK1rankeduniqueresult.csv")
rownames(PINK1) <- PINK1$Gene.Symbol
PINK1 <- PINK1[,204:215]
PINK1 <- subset(PINK1, rownames(PINK1) %in% DEG_PPI)


# Cor.test Method #

library(tictoc)
library(gdata)

##For loop for generating regression values and p values
CorExprMat <- t(PINK1)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
p.value <- matrix(0, ncol(test), ncol(test))

Sys.time()
tic()
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(test)
toc()
Sys.time()
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
write.csv(results, "PINK1_coexpression.csv")


## filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")

ATP13A2 <- read.csv("ATP13A2_coexpression.csv")
ATP13A2$Gene1 <- as.character(lapply(strsplit(as.character(ATP13A2$X), "\\:"), "[", 2))
ATP13A2$Gene2 <- as.character(lapply(strsplit(as.character(ATP13A2$X), "\\:"), "[", 1))
ATP13A2 <- ATP13A2[,c(5,6,1,2,3,4)]

PARKIN <- read.csv("PARKIN_coexpression.csv")
PARKIN$Gene1 <- as.character(lapply(strsplit(as.character(PARKIN$X), "\\:"), "[", 2))
PARKIN$Gene2 <- as.character(lapply(strsplit(as.character(PARKIN$X), "\\:"), "[", 1))
PARKIN <- PARKIN[,c(5,6,1,2,3,4)]

PINK1 <- read.csv("PINK1_coexpression.csv")
PINK1$Gene1 <- as.character(lapply(strsplit(as.character(PINK1$X), "\\:"), "[", 2))
PINK1$Gene2 <- as.character(lapply(strsplit(as.character(PINK1$X), "\\:"), "[", 1))
PINK1 <- PINK1[,c(5,6,1,2,3,4)]

hist(ATP13A2$reg.mat)
hist(PARKIN$reg.mat)
hist(PINK1$reg.mat)

thresh <- 0.5
### Filter by r value
ATP13A2_cor.5 <- ATP13A2[ATP13A2$reg.mat > thresh | ATP13A2$reg.mat < -thresh,]
PARKIN_cor.5 <- PARKIN[PARKIN$reg.mat > thresh | PARKIN$reg.mat < -thresh,]
PINK1_cor.5 <- PINK1[PINK1$reg.mat > thresh | PINK1$reg.mat < -thresh,]



ATP13A2_cor.5$X2 <- paste(ATP13A2_cor.5$Gene1,":",ATP13A2_cor.5$Gene2, sep = "")
PARKIN_cor.5$X2 <- paste(PARKIN_cor.5$Gene1,":",PARKIN_cor.5$Gene2, sep = "")
PINK1_cor.5$X2 <- paste(PINK1_cor.5$Gene1,":",PINK1_cor.5$Gene2, sep = "")


#Merge first two datasets
corresult1 <- merge(ATP13A2_cor.5, PARKIN_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(ATP13A2_cor.5, PARKIN_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "ATP13A2", "PARKIN")

#Merge with 3rd Dataset
corresult3 <- merge(result, PINK1_cor.5, by.x = "Correlation", by.y = "X")
corresult4 <- merge(result, PINK1_cor.5, by.x = "Correlation", by.y = "X2")

corresult3 <- corresult3[, c(1:5, 10)]
corresult4 <- corresult4[, c(1:5, 11)]

result <- rbind(corresult3, corresult4)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "ATP13A2", "PARKIN", "PINK1")



CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:6]

CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "famblood_0point5_2.csv", quote = F, row.names = F)


pat <- read.csv("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PD_0point5_wDUM.csv")
con <- read.csv("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/Control_0point5_2.csv")

patcon <- merge(pat, con, by = "Gene")
patcon <- patcon[,c(1,9,20)]
patcon$Gene1 <- as.character(lapply(strsplit(as.character(patcon$Gene), "\\:"), "[", 2))
patcon$Gene2 <- as.character(lapply(strsplit(as.character(patcon$Gene), "\\:"), "[", 1))
write.csv(patcon, "overlap_remove.csv", row.names = F, quote = F)

####################################
#### ALS Network ####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PD_NuggetGenes_conoverlapremoved.txt")


setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/")

C9 <- read.csv("C9uniquegene_samples.csv", row.names = 1)
C9 <- C9[,4:11]
C9 <- subset(C9, rownames(C9) %in% DEG_PPI)

sals <- read.csv("sALSuniquegene_samples.csv", row.names = 1)
sals <- sals[,4:10]
sals <- subset(sals, rownames(sals) %in% DEG_PPI)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_results_keepfiltering.csv")
rownames(pet) <- pet$hgnc_symbol
pet <- pet[,19:35]
pet <- subset(pet, rownames(pet) %in% DEG_PPI)

rav <- read.csv("RAV_results_keepfiltering.csv")
rownames(rav) <- rav$hgnc_symbol
rav <- rav[,17:29]
rav <- subset(rav, rownames(rav) %in% DEG_PPI)


# Cor.test Method #

library(tictoc)
library(gdata)

##For loop for generating regression values and p values
CorExprMat <- t(rav)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
p.value <- matrix(0, ncol(test), ncol(test))

Sys.time()
tic()
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(test)
toc()
Sys.time()
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
write.csv(results, "rav_coexpression.csv")


## filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")

C9 <- read.csv("C9_coexpression.csv")
C9$Gene1 <- as.character(lapply(strsplit(as.character(C9$X), "\\:"), "[", 2))
C9$Gene2 <- as.character(lapply(strsplit(as.character(C9$X), "\\:"), "[", 1))
C9 <- C9[,c(5,6,1,2,3,4)]

sals <- read.csv("sals_coexpression.csv")
sals$Gene1 <- as.character(lapply(strsplit(as.character(sals$X), "\\:"), "[", 2))
sals$Gene2 <- as.character(lapply(strsplit(as.character(sals$X), "\\:"), "[", 1))
sals <- sals[,c(5,6,1,2,3,4)]

pet <- read.csv("pet_coexpression.csv")
pet$Gene1 <- as.character(lapply(strsplit(as.character(pet$X), "\\:"), "[", 2))
pet$Gene2 <- as.character(lapply(strsplit(as.character(pet$X), "\\:"), "[", 1))
pet <- pet[,c(5,6,1,2,3,4)]

rav <- read.csv("rav_coexpression.csv")
rav$Gene1 <- as.character(lapply(strsplit(as.character(rav$X), "\\:"), "[", 2))
rav$Gene2 <- as.character(lapply(strsplit(as.character(rav$X), "\\:"), "[", 1))
rav <- rav[,c(5,6,1,2,3,4)]

hist(C9$reg.mat)
hist(sals$reg.mat)
hist(pet$reg.mat)
hist(rav$reg.mat)

thresh <- 0.5
### Filter by r value
C9_cor.5 <- C9[C9$reg.mat > thresh | C9$reg.mat < -thresh,]
sals_cor.5 <- sals[sals$reg.mat > thresh | sals$reg.mat < -thresh,]
pet_cor.5 <- pet[pet$reg.mat > thresh | pet$reg.mat < -thresh,]
rav_cor.5 <- rav[rav$reg.mat > thresh | rav$reg.mat < -thresh,]



C9_cor.5$X2 <- paste(C9_cor.5$Gene1,":",C9_cor.5$Gene2, sep = "")
sals_cor.5$X2 <- paste(sals_cor.5$Gene1,":",sals_cor.5$Gene2, sep = "")
pet_cor.5$X2 <- paste(pet_cor.5$Gene1,":",pet_cor.5$Gene2, sep = "")
rav_cor.5$X2 <- paste(rav_cor.5$Gene1,":",rav_cor.5$Gene2, sep = "")


#Merge first two datasets
corresult1 <- merge(C9_cor.5, sals_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(C9_cor.5, sals_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sals")

#Merge with 3rd Dataset
corresult3 <- merge(result, pet_cor.5, by.x = "Correlation", by.y = "X")
corresult4 <- merge(result, pet_cor.5, by.x = "Correlation", by.y = "X2")

corresult3 <- corresult3[, c(1:5, 10)]
corresult4 <- corresult4[, c(1:5, 11)]

result <- rbind(corresult3, corresult4)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sals", "pet")

#Merge with 4th Dataset
corresult5 <- merge(result, rav_cor.5, by.x = "Correlation", by.y = "X")
corresult6 <- merge(result, rav_cor.5, by.x = "Correlation", by.y = "X2")

corresult5 <- corresult5[, c(1:6, 11)]
corresult6 <- corresult6[, c(1:6, 12)]

result <- rbind(corresult5, corresult6)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "C9", "sals", "pet", "rav")

CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:7]

CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "ALS_0point5.csv", quote = F, row.names = F)



setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")



####################################
#### Sporadic Blood Network ####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PD_NuggetGenes_conoverlapremoved.txt")

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
AMA <- read.csv("AMAfilteredresult.csv")
rownames(AMA) <- AMA$Gene.Symbol
AMA <- AMA[,204:372]
AMA <- subset(AMA, rownames(AMA) %in% DEG_PPI)

RON <- read.csv("RONfilteredresult.csv")
rownames(RON) <- RON$Gene.Symbol
RON <- RON[,39:78]
RON <- subset(RON, rownames(RON) %in% DEG_PPI)



# Cor.test Method #

library(tictoc)
library(gdata)

##For loop for generating regression values and p values
CorExprMat <- t(RON)
test <- CorExprMat
reg <- matrix(0, ncol(test), ncol(test))
p.value <- matrix(0, ncol(test), ncol(test))

Sys.time()
tic()
for (i in 1:ncol(test)){
  for (j in 1:ncol(test)){
    reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(test)
toc()
Sys.time()
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

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
write.csv(results, "RON_coexpression.csv")


## filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")

AMA <- read.csv("AMA_coexpression.csv")
AMA$Gene1 <- as.character(lapply(strsplit(as.character(AMA$X), "\\:"), "[", 2))
AMA$Gene2 <- as.character(lapply(strsplit(as.character(AMA$X), "\\:"), "[", 1))
AMA <- AMA[,c(5,6,1,2,3,4)]

RON <- read.csv("RON_coexpression.csv")
RON$Gene1 <- as.character(lapply(strsplit(as.character(RON$X), "\\:"), "[", 2))
RON$Gene2 <- as.character(lapply(strsplit(as.character(RON$X), "\\:"), "[", 1))
RON <- RON[,c(5,6,1,2,3,4)]


hist(AMA$reg.mat)
hist(RON$reg.mat)


thresh <- 0.5
### Filter by r value
AMA_cor.5 <- AMA[AMA$reg.mat > thresh | AMA$reg.mat < -thresh,]
RON_cor.5 <- RON[RON$reg.mat > thresh | RON$reg.mat < -thresh,]


AMA_cor.5$X2 <- paste(AMA_cor.5$Gene1,":",AMA_cor.5$Gene2, sep = "")
RON_cor.5$X2 <- paste(RON_cor.5$Gene1,":",RON_cor.5$Gene2, sep = "")


#Merge first two datasets
corresult1 <- merge(AMA_cor.5, RON_cor.5, by.x = "X", by.y = "X")
corresult2 <- merge(AMA_cor.5, RON_cor.5, by.x = "X2", by.y = "X")

#Only take relevent columns
corresult1 <- corresult1[, c(1:3, 6, 12)]
corresult2 <- corresult2[, c(1:3, 7, 12)]
colnames(corresult2)[1] <- "X"

#bind one under the other
result <- rbind(corresult1, corresult2)
colnames(result) <- c("Correlation", "Gene1", "Gene2", "AMA", "RON")

CommonGroup <- result
CommonGroup <- CommonGroup[!duplicated(CommonGroup[,1]),]
rownames(CommonGroup) <- CommonGroup$Correlation
CommonGroup <- CommonGroup[,4:5]

CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ControlCoexpression/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "spblood_0point5.csv", quote = F, row.names = F)
