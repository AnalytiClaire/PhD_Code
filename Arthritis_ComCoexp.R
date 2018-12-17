#### Arthritis DEG PPI Correlation ####

#Read in network nodes

#######
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PPI/ArthritisPPI.txt")

setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
#Extract PPI network genes from each dataset#####
BER.OA <- read.csv("BER_OA_filteredresult.csv")
rownames(BER.OA) <- BER.OA$Gene.Symbol
BER.OA <- BER.OA[,30:39]
BER.OA <- subset(BER.OA, rownames(BER.OA) %in% DEG_PPI)

BER.RA <- read.csv("BER_RA_filteredresult.csv")
rownames(BER.RA) <- BER.RA$Gene.Symbol
BER.RA <- BER.RA[,30:39]
BER.RA <- subset(BER.RA, rownames(BER.RA) %in% DEG_PPI)

JEN.OA <- read.csv("JEN_OA_filteredresult.csv")
rownames(JEN.OA) <- JEN.OA$Gene.Symbol
JEN.OA <- JEN.OA[,30:39]
JEN.OA <- subset(JEN.OA, rownames(JEN.OA) %in% DEG_PPI)

JEN.RA <- read.csv("JEN_RA_filteredresult.csv")
rownames(JEN.RA) <- JEN.RA$Gene.Symbol
JEN.RA <- JEN.RA[,30:42]
JEN.RA <- subset(JEN.RA, rownames(JEN.RA) %in% DEG_PPI)

BRO <- read.csv("BRO_RA_filteredresult.csv")
rownames(BRO) <- BRO$Gene.Symbol
BRO <- BRO[,28:43]
BRO <- subset(BRO, rownames(BRO) %in% DEG_PPI)

VRI <- read.csv("VRI_OA_filteredresult.csv")
rownames(VRI) <- VRI$Gene.Symbol
VRI <- VRI[,28:37]
VRI <- subset(VRI, rownames(VRI) %in% DEG_PPI)

PAD <- read.csv("PAD_uniqueresult.csv")
rownames(PAD) <- PAD$hgnc_symbol
PAD <- PAD[,17:25]
PAD <- subset(PAD, rownames(PAD) %in% DEG_PPI)

WAL.OA <- read.csv("WAL_OA_UniqueGene_DESeq2.csv")
rownames(WAL.OA) <- WAL.OA$Row.names
WAL.OA <- WAL.OA[,37:58]
WAL.OA <- subset(WAL.OA, rownames(WAL.OA) %in% DEG_PPI)

WAL.RA <- read.csv("WAL_RA_UniqueGene_DESeq2.csv")
rownames(WAL.RA) <- WAL.RA$Row.names
WAL.RA <- WAL.RA[,37:183]
WAL.RA <- subset(WAL.RA, rownames(WAL.RA) %in% DEG_PPI)




#Find the gene names that all datasets have in common
DEG_com <- Reduce(intersect, list(rownames(BER.OA), rownames(BER.RA), 
                                  rownames(BRO),rownames(JEN.OA),rownames(JEN.RA),
                                  rownames(PAD),rownames(VRI),rownames(WAL.OA),
                                  rownames(WAL.RA)))

#Subset each dataset with these common names so they are all the same size
BER.OA <- subset(BER.OA, rownames(BER.OA) %in% DEG_com)
BER.RA <- subset(BER.RA, rownames(BER.RA) %in% DEG_com)
BRO <- subset(BRO, rownames(BRO) %in% DEG_com)
JEN.OA <- subset(JEN.OA, rownames(JEN.OA) %in% DEG_com)
JEN.RA <- subset(JEN.RA, rownames(JEN.RA) %in% DEG_com)
PAD <- subset(PAD, rownames(PAD) %in% DEG_com)
VRI <- subset(VRI, rownames(VRI) %in% DEG_com)
WAL.OA <- subset(WAL.OA, rownames(WAL.OA) %in% DEG_com)
WAL.RA <- subset(WAL.RA, rownames(WAL.RA) %in% DEG_com)

setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/ComCoexpr/")

write.csv(BER.OA, "BER.OA_DEGfilt.csv")
write.csv(BER.RA, "BER.RA_DEGfilt.csv")
write.csv(BRO, "BRO_DEGfilt.csv")
write.csv(JEN.OA, "JEN.OA_DEGfilt.csv")
write.csv(JEN.RA, "JEN.RA_DEGfilt.csv")
write.csv(PAD, "PAD_DEGfilt.csv")
write.csv(VRI, "VRI_DEGfilt.csv")
write.csv(WAL.OA, "WAL.OA_DEGfilt.csv")
write.csv(WAL.RA, "WAL.RA_DEGfilt.csv")

#Run Correlation analysis on Sharc


#### filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/ComCoexpr/")

BER.OA <- read.csv("BER.OAcorresult.csv")
BER.OA$Gene1 <- as.character(lapply(strsplit(as.character(BER.OA$X), "\\:"), "[", 2))
BER.OA$Gene2 <- as.character(lapply(strsplit(as.character(BER.OA$X), "\\:"), "[", 1))
BER.OA <- BER.OA[,c(5,6,1,2,3,4)]

BER.RA <- read.csv("BER.RAcorresult.csv")
BER.RA$Gene1 <- as.character(lapply(strsplit(as.character(BER.RA$X), "\\:"), "[", 2))
BER.RA$Gene2 <- as.character(lapply(strsplit(as.character(BER.RA$X), "\\:"), "[", 1))
BER.RA <- BER.RA[,c(5,6,1,2,3,4)]

BRO <- read.csv("BROcorresult.csv")
BRO$Gene1 <- as.character(lapply(strsplit(as.character(BRO$X), "\\:"), "[", 2))
BRO$Gene2 <- as.character(lapply(strsplit(as.character(BRO$X), "\\:"), "[", 1))
BRO <- BRO[,c(5,6,1,2,3,4)]

JEN.OA <- read.csv("JEN.OAcorresult.csv")
JEN.OA$Gene1 <- as.character(lapply(strsplit(as.character(JEN.OA$X), "\\:"), "[", 2))
JEN.OA$Gene2 <- as.character(lapply(strsplit(as.character(JEN.OA$X), "\\:"), "[", 1))
JEN.OA <- JEN.OA[,c(5,6,1,2,3,4)]

JEN.RA <- read.csv("JEN.RAcorresult.csv")
JEN.RA$Gene1 <- as.character(lapply(strsplit(as.character(JEN.RA$X), "\\:"), "[", 2))
JEN.RA$Gene2 <- as.character(lapply(strsplit(as.character(JEN.RA$X), "\\:"), "[", 1))
JEN.RA <- JEN.RA[,c(5,6,1,2,3,4)]

PAD <- read.csv("PADcorresult.csv")
PAD$Gene1 <- as.character(lapply(strsplit(as.character(PAD$X), "\\:"), "[", 2))
PAD$Gene2 <- as.character(lapply(strsplit(as.character(PAD$X), "\\:"), "[", 1))
PAD <- PAD[,c(5,6,1,2,3,4)]
PAD$X <- paste(PAD$Gene1,":",PAD$Gene2, sep = "")

VRI <- read.csv("VRIcorresult.csv")
VRI$Gene1 <- as.character(lapply(strsplit(as.character(VRI$X), "\\:"), "[", 2))
VRI$Gene2 <- as.character(lapply(strsplit(as.character(VRI$X), "\\:"), "[", 1))
VRI <- VRI[,c(5,6,1,2,3,4)]
VRI$X <- paste(VRI$Gene1,":",VRI$Gene2, sep = "")


WAL.OA <- read.csv("WAL.OAcorresult.csv")
WAL.OA$Gene1 <- as.character(lapply(strsplit(as.character(WAL.OA$X), "\\:"), "[", 2))
WAL.OA$Gene2 <- as.character(lapply(strsplit(as.character(WAL.OA$X), "\\:"), "[", 1))
WAL.OA <- WAL.OA[,c(5,6,1,2,3,4)]
WAL.OA$X <- paste(WAL.OA$Gene2,":",WAL.OA$Gene1, sep = "")

WAL.RA<- read.csv("WAL.RAcorresult.csv")
WAL.RA$Gene1 <- as.character(lapply(strsplit(as.character(WAL.RA$X), "\\:"), "[", 2))
WAL.RA$Gene2 <- as.character(lapply(strsplit(as.character(WAL.RA$X), "\\:"), "[", 1))
WAL.RA <- WAL.RA[,c(5,6,1,2,3,4)]


thresh <- 0.3
### Filter by r value
BER.OA_cor.5 <- BER.OA[BER.OA$reg.mat > thresh | BER.OA$reg.mat < -thresh,]
BER.RA_cor.5 <- BER.RA[BER.RA$reg.mat > thresh | BER.RA$reg.mat < -thresh,]
BRO_cor.5    <- BRO[BRO$reg.mat > thresh | BRO$reg.mat < -thresh,]
JEN.OA_cor.5 <- JEN.OA[JEN.OA$reg.mat > thresh | JEN.OA$reg.mat < -thresh,]
JEN.RA_cor.5 <- JEN.RA[JEN.RA$reg.mat > thresh | JEN.RA$reg.mat < -thresh,]
# PAD_cor.5    <- PAD[PAD$reg.mat > thresh | PAD$reg.mat < -thresh,]
VRI_cor.5    <- VRI[VRI$reg.mat > thresh | VRI$reg.mat < -thresh,]
# WAL.OA_cor.5 <- WAL.OA[WAL.OA$reg.mat > thresh | WAL.OA$reg.mat < -thresh,]
WAL.RA_cor.5 <- WAL.RA[WAL.RA$reg.mat > thresh | WAL.RA$reg.mat < -thresh,]


### Find matches 

Commonedge <- Reduce(intersect, list(BER.OA_cor.5$X, BER.RA_cor.5$X, BRO_cor.5$X, 
                                     JEN.OA_cor.5$X, JEN.RA_cor.5$X,
                                     VRI_cor.5$X, WAL.RA_cor.5$X))


#Subset each dataset with these common names so they are all the same size
BER.OA_CE <- subset(BER.OA_cor.5, BER.OA_cor.5$X %in% Commonedge)
BER.OA_CE <- BER.OA_CE[order(BER.OA_CE$X),]
BER.RA_CE <- subset(BER.RA_cor.5, BER.RA_cor.5$X %in% Commonedge)
BER.RA_CE <- BER.RA_CE[order(BER.RA_CE$X),]
BRO_CE <- subset(BRO_cor.5, BRO_cor.5$X %in% Commonedge)
BRO_CE <- BRO_CE[order(BRO_CE$X),]
JEN.OA_CE <- subset(JEN.OA_cor.5, JEN.OA_cor.5$X %in% Commonedge)
JEN.OA_CE <- JEN.OA_CE[order(JEN.OA_CE$X),]
JEN.RA_CE <- subset(JEN.RA_cor.5, JEN.RA_cor.5$X %in% Commonedge)
JEN.RA_CE <- JEN.RA_CE[order(JEN.RA_CE$X),]
# PAD_CE <- subset(PAD_cor.5, PAD_cor.5$X %in% Commonedge)
# PAD_CE <- PAD_CE[order(PAD_CE$X),]
VRI_CE <- subset(VRI_cor.5, VRI_cor.5$X %in% Commonedge)
VRI_CE <- VRI_CE[order(VRI_CE$X),]
# WAL.OA_CE <- subset(WAL.OA_cor.5, WAL.OA_cor.5$X %in% Commonedge)
# WAL.OA_CE <- WAL.OA_CE[order(WAL.OA_CE$X),]
WAL.RA_CE <- subset(WAL.RA_cor.5, WAL.RA_cor.5$X %in% Commonedge)
WAL.RA_CE <- WAL.RA_CE[order(WAL.RA_CE$X),]



CommonGroup <- data.frame(row.names = BER.OA_CE$X,
                          BER.OA = BER.OA_CE$reg.mat,
                          BER.RA = BER.RA_CE$reg.mat,
                          BRO = BRO_CE$reg.mat,
                          JEN.OA = JEN.OA_CE$reg.mat,
                          JEN.RA = JEN.RA_CE$reg.mat,
                          VRI = VRI_CE$reg.mat,
                          WAL.RA = WAL.RA_CE$reg.mat)


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/ComCoexpr/")
CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "RAOA_0point3_NoWAL.OAPAD.csv", quote = F, row.names = F)

nodetable <- read.csv("RAOA_0point3_NoWAL.OAPAD.nuggetnode.csv")
nuggenes <- as.character(nodetable$name)
RAGenes <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/RAGenes.txt")
OAGenes <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/OAGenes.txt")
celltype <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv")
DEG <- readLines("../GeneExpression/allDEGs_SJO_SLE_ALO.txt")
Okada <- readLines("../../Okada_MetaGWAS_RA.txt")

# nalls <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/NallsPDGWAS.txt")
# targetvalPD <- c("SIRT2", "DNM1L", "STMN1", "DNAJA3", "TBK1", "RTCA", "ANXA7", "DNAJC12", "RTN4", "ADSL", "MDH1","ATP6V1G2", 
#                  "YWHAZ", "ETS1")

nodetable$RAGenes <- nodetable$shared.name %in% RAGenes
nodetable$OAGenes <- nodetable$shared.name %in% OAGenes
nodetable$DEG <- nodetable$name %in% DEG
nodetable$OkadaGWAS <- nodetable$name %in% Okada
# nodetable$targetvalPD <- nodetable$shared.name %in% targetvalPD
# nodetable$targetvalPD <- nodetable$shared.name %in% nalls

nodetable_celltype <- merge(celltype, nodetable, by.x = "Gene.symbol",  by.y = "shared.name", all = T)
nodetable_celltype <- subset(nodetable_celltype, !(nodetable_celltype$SUID == "NA"))


# nodetable_celltype$targetvalPD <- nodetable_celltype$Gene.symbol %in% targetvalPD
write.csv(nodetable_celltype, "ModifiedRAOANodetable.csv", row.names = F)

### Control network ####
#### Familial Blood ####
#Read in network nodes
DEG_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/0point3nuggetgenes.txt")

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

DIJ <- read.csv("DIJfilteredresult.csv")
rownames(DIJ) <- DIJ$Gene.Symbol
DIJ <- DIJ[,21:28]
DIJ <- subset(DIJ, rownames(DIJ) %in% DEG_PPI)

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

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/ControlCoexpression/")

write.csv(DIJ, "DIJ_DEGfilt.csv")
write.csv(FFR, "FFR_DEGfilt.csv")
write.csv(MID1, "MID1_DEGfilt.csv")
write.csv(MID2, "MID2_DEGfilt.csv")
write.csv(MID3, "MID3_DEGfilt.csv")
write.csv(MID4, "MID4_DEGfilt.csv")
write.csv(MOR.SN, "MOR.SN_DEGfilt.csv")


#### filter correlations #### BOT/BOT2 not included due to less than 4 samples
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/ControlCoexpression/")

DIJ <- read.csv("DIJcorresult.csv")
DIJ$Gene1 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 2))
DIJ$Gene2 <- as.character(lapply(strsplit(as.character(DIJ$X), "\\:"), "[", 1))
DIJ <- DIJ[,c(5,6,1,2,3,4)]

FFR <- read.csv("FFRcorresult.csv")
FFR$Gene1 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 2))
FFR$Gene2 <- as.character(lapply(strsplit(as.character(FFR$X), "\\:"), "[", 1))
FFR <- FFR[,c(5,6,1,2,3,4)]
FFR$X <- paste(FFR$Gene1,":",FFR$Gene2, sep = "")

MID1 <- read.csv("MID1corresult.csv")
MID1$Gene1 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 2))
MID1$Gene2 <- as.character(lapply(strsplit(as.character(MID1$X), "\\:"), "[", 1))
MID1 <- MID1[,c(5,6,1,2,3,4)]

MID2 <- read.csv("MID2corresult.csv")
MID2$Gene1 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 2))
MID2$Gene2 <- as.character(lapply(strsplit(as.character(MID2$X), "\\:"), "[", 1))
MID2 <- MID2[,c(5,6,1,2,3,4)]

MID3 <- read.csv("MID3corresult.csv")
MID3$Gene1 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 2))
MID3$Gene2 <- as.character(lapply(strsplit(as.character(MID3$X), "\\:"), "[", 1))
MID3 <- MID3[,c(5,6,1,2,3,4)]
MID3$X <- paste(MID3$Gene1,":",MID3$Gene2, sep = "")

MID4 <- read.csv("MID4corresult.csv")
MID4$Gene1 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 2))
MID4$Gene2 <- as.character(lapply(strsplit(as.character(MID4$X), "\\:"), "[", 1))
MID4 <- MID4[,c(5,6,1,2,3,4)]

MOR.SN<- read.csv("MOR.SNcorresult.csv")
MOR.SN$Gene1 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 2))
MOR.SN$Gene2 <- as.character(lapply(strsplit(as.character(MOR.SN$X), "\\:"), "[", 1))
MOR.SN <- MOR.SN[,c(5,6,1,2,3,4)]
MOR.SN$X <- paste(MOR.SN$Gene1,":",MOR.SN$Gene2, sep = "")

thresh <- 0.3
### Filter by r value
DIJ_cor.5 <- DIJ[DIJ$reg.mat > thresh | DIJ$reg.mat < -thresh,]
FFR_cor.5 <- FFR[FFR$reg.mat > thresh | FFR$reg.mat < -thresh,]
MID1_cor.5 <- MID1[MID1$reg.mat > thresh | MID1$reg.mat < -thresh,]
MID2_cor.5 <- MID2[MID2$reg.mat > thresh | MID2$reg.mat < -thresh,]
MID3_cor.5 <- MID3[MID3$reg.mat > thresh | MID3$reg.mat < -thresh,]
MID4_cor.5 <- MID4[MID4$reg.mat > thresh | MID4$reg.mat < -thresh,]
MOR.SN_cor.5 <- MOR.SN[MOR.SN$reg.mat > thresh | MOR.SN$reg.mat < -thresh,]


### Find matches 

Commonedge <- Reduce(intersect, list(DIJ_cor.5$X,FFR_cor.5$X, MID1_cor.5$X, 
                                     MID2_cor.5$X,MID3_cor.5$X, MID4_cor.5$X,
                                     MOR.SN_cor.5$X))


#Subset each dataset with these common names so they are all the same size
DIJ_CE <- subset(DIJ_cor.5, DIJ_cor.5$X %in% Commonedge)
DIJ_CE <- DIJ_CE[order(DIJ_CE$X),]
FFR_CE <- subset(FFR_cor.5, FFR_cor.5$X %in% Commonedge)
FFR_CE <- FFR_CE[order(FFR_CE$X),]
MID1_CE <- subset(MID1_cor.5, MID1_cor.5$X %in% Commonedge)
MID1_CE <- MID1_CE[order(MID1_CE$X),]
MID2_CE <- subset(MID2_cor.5, MID2_cor.5$X %in% Commonedge)
MID2_CE <- MID2_CE[order(MID2_CE$X),]
MID3_CE <- subset(MID3_cor.5, MID3_cor.5$X %in% Commonedge)
MID3_CE <- MID3_CE[order(MID3_CE$X),]
MID4_CE <- subset(MID4_cor.5, MID4_cor.5$X %in% Commonedge)
MID4_CE <- MID4_CE[order(MID4_CE$X),]
MOR.SN_CE <- subset(MOR.SN_cor.5, MOR.SN_cor.5$X %in% Commonedge)
MOR.SN_CE <- MOR.SN_CE[order(MOR.SN_CE$X),]



CommonGroup <- data.frame(row.names = DIJ_CE$X,
                          DIJ = DIJ_CE$reg.mat,
                          FFR = FFR_CE$reg.mat,
                          MID1 = MID1_CE$reg.mat,
                          MID2 = MID2_CE$reg.mat, 
                          MID3 = MID3_CE$reg.mat, 
                          MID4 = MID4_CE$reg.mat, 
                          MOR.SN = MOR.SN_CE$reg.mat)


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
CG_samedir$Gene <- rownames(CG_samedir)
CG_samedir$Gene1 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 2))
CG_samedir$Gene2 <- as.character(lapply(strsplit(as.character(CG_samedir$Gene), "\\:"), "[", 1))

write.csv(CG_samedir, "PD_0point3.csv", quote = F, row.names = F)
print(CG_conserved_up)
