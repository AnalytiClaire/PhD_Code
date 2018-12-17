#### Parkinson's DEG PPI Correlation ####

#Read in network nodes

#### ALP####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/ALSPARK2LRRK2/ALP_PPIGenes.txt")

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
#Extract PPI network genes from each dataset#####
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

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/ALSPARK2LRRK2/")

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


#### filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/ALSPARK2LRRK2/ComCoexp/")

DIJ <- read.csv("DIJcorresult.csv")
DIJ$Gene1 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 2))
DIJ$Gene2 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 1))
DIJ <- DIJ[,c(5,6,1,2,3,4)]

DUM <- read.csv("DUMcorresult.csv")
DUM$Gene1 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 2))
DUM$Gene2 <- as.character(lapply(strsplit(as.character(DUM$X), "\\:"), "[", 1))
DUM <- DUM[,c(5,6,1,2,3,4)]
DUM$X <- paste(DUM$Gene1,":",DUM$Gene2, sep = "")

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
MID4$X <- paste(MID4$Gene1,":",MID4$Gene2, sep = "")

MOR.FC<- read.csv("MOR.FCcorresult.csv")
MOR.FC$Gene1 <- as.character(lapply(strsplit(as.character(MOR.FC$X), "\\:"), "[", 2))
MOR.FC$Gene2 <- as.character(lapply(strsplit(as.character(MOR.FC$X), "\\:"), "[", 1))
MOR.FC <- MOR.FC[,c(5,6,1,2,3,4)]

MOR.SN<- read.csv("MOR.SNcorresult.csv")
MOR.SN$Gene1 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 2))
MOR.SN$Gene2 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 1))
MOR.SN <- MOR.SN[,c(5,6,1,2,3,4)]

thresh <- 0.1
### Filter by r value
DIJ_cor.5 <- DIJ[DIJ$reg.mat > thresh | DIJ$reg.mat < -thresh,]
DUM_cor.5 <- DUM[DUM$reg.mat > thresh | DUM$reg.mat < -thresh,]
FFR_cor.5 <- FFR[FFR$reg.mat > thresh | FFR$reg.mat < -thresh,]
LEW_cor.5 <- LEW[LEW$reg.mat > thresh | LEW$reg.mat < -thresh,]
MID1_cor.5 <- MID1[MID1$reg.mat > thresh | MID1$reg.mat < -thresh,]
MID2_cor.5 <- MID2[MID2$reg.mat > thresh | MID2$reg.mat < -thresh,]
MID3_cor.5 <- MID3[MID3$reg.mat > thresh | MID3$reg.mat < -thresh,]
MID4_cor.5 <- MID4[MID4$reg.mat > thresh | MID4$reg.mat < -thresh,]
MOR.FC_cor.5 <- MOR.FC[MOR.FC$reg.mat > thresh | MOR.FC$reg.mat < -thresh,]
MOR.SN_cor.5 <- MOR.SN[MOR.SN$reg.mat > thresh | MOR.SN$reg.mat < -thresh,]


### Find matches 

Commonedge <- Reduce(intersect, list(DIJ_cor.5$X, FFR_cor.5$X,
                                     LEW_cor.5$X, MID1_cor.5$X, MID2_cor.5$X,
                                     MID3_cor.5$X, MID4_cor.5$X,
                                     MOR.SN_cor.5$X))


#Subset each dataset with these common names so they are all the same size
DIJ_CE <- subset(DIJ_cor.5, DIJ_cor.5$X %in% Commonedge)
DIJ_CE <- DIJ_CE[order(DIJ_CE$X),]
# DUM_CE <- subset(DUM_cor.5, DUM_cor.5$X %in% Commonedge)
# DUM_CE <- DUM_CE[order(DUM_CE$X),]
FFR_CE <- subset(FFR_cor.5, FFR_cor.5$X %in% Commonedge)
FFR_CE <- FFR_CE[order(FFR_CE$X),]
LEW_CE <- subset(LEW_cor.5, LEW_cor.5$X %in% Commonedge)
LEW_CE <- LEW_CE[order(LEW_CE$X),]
MID1_CE <- subset(MID1_cor.5, MID1_cor.5$X %in% Commonedge)
MID1_CE <- MID1_CE[order(MID1_CE$X),]
MID2_CE <- subset(MID2_cor.5, MID2_cor.5$X %in% Commonedge)
MID2_CE <- MID2_CE[order(MID2_CE$X),]
MID3_CE <- subset(MID3_cor.5, MID3_cor.5$X %in% Commonedge)
MID3_CE <- MID3_CE[order(MID3_CE$X),]
MID4_CE <- subset(MID4_cor.5, MID4_cor.5$X %in% Commonedge)
MID4_CE <- MID4_CE[order(MID4_CE$X),]
# MOR.FC_CE <- subset(MOR.FC_cor.5, MOR.FC_cor.5$X %in% Commonedge)
# MOR.FC_CE <- MOR.FC_CE[order(MOR.FC_CE$X),]
MOR.SN_CE <- subset(MOR.SN_cor.5, MOR.SN_cor.5$X %in% Commonedge)
MOR.SN_CE <- MOR.SN_CE[order(MOR.SN_CE$X),]



CommonGroup <- data.frame(row.names = DIJ_CE$X,
                          DIJ = DIJ_CE$reg.mat,
                          FFR = FFR_CE$reg.mat,
                          LEW = LEW_CE$reg.mat,
                          MID1 = MID1_CE$reg.mat,
                          MID2 = MID2_CE$reg.mat, 
                          MID3 = MID3_CE$reg.mat, 
                          MID4 = MID4_CE$reg.mat, 
                          MOR.SN = MOR.SN_CE$reg.mat)


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/ALSPARK2LRRK2/ComCoexp/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))



write.csv(CG_samedir, "RMETHOD_samedir_mean.csv", quote = F, row.names = F)
