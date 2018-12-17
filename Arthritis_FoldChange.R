setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")

BER.OA <- read.csv("BER_OA_filteredresult.csv")
BER.OA <- BER.OA[order(BER.OA$P.Value),]
BER.RA <- read.csv("BER_RA_filteredresult.csv")
BER.RA <- BER.RA[order(BER.RA$P.Value),]
JEN.OA <- read.csv("JEN_OA_filteredresult.csv")
JEN.OA <- JEN.OA[order(JEN.OA$P.Value),]
JEN.RA <- read.csv("JEN_RA_filteredresult.csv")
JEN.RA <- JEN.RA[order(JEN.RA$P.Value),]
BRO <- read.csv("BRO_RA_filteredresult.csv")
BRO <- BRO[order(BRO$P.Value),]
VRI <- read.csv("VRI_OA_filteredresult.csv")
VRI <- VRI[order(VRI$P.Value),]
PAD <- read.csv("PAD_uniqueresult.csv")
PAD <- PAD[order(PAD$pvalue),]
WAL.OA <- read.csv("WAL_OA_UniqueGene_DESeq2.csv")
WAL.OA <- WAL.OA[order(WAL.OA$pvalue),]
WAL.RA <- read.csv("WAL_RA_UniqueGene_DESeq2.csv")
WAL.RA <- WAL.RA[order(WAL.RA$pvalue),]




thresh <- 1

upBER.OA <- subset(BER.OA, BER.OA$Fold.Change >= thresh)
upBER.OAgene <- upBER.OA$Gene.Symbol

upBER.RA <- subset(BER.RA, BER.RA$Fold.Change >= thresh)
upBER.RAgene <- upBER.RA$Gene.Symbol

upJEN.OA <- subset(JEN.OA, JEN.OA$Fold.Change >= thresh)
upJEN.OAgene <- upJEN.OA$Gene.Symbol

upJEN.RA <- subset(JEN.RA, JEN.RA$Fold.Change >= thresh)
upJEN.RAgene <- upJEN.RA$Gene.Symbol

upBRO <- subset(BRO, BRO$Fold.Change >= thresh)
upBROgene <- upBRO$Gene.Symbol

upVRI <- subset(VRI, VRI$Fold.Change >= thresh)
upVRIgene <- upVRI$Gene.Symbol

upPAD <- subset(PAD, PAD$Fold.Change >= thresh)
upPADgene <- upPAD$hgnc_symbol

upWAL.OA <- subset(WAL.OA, WAL.OA$Fold.Change >= thresh)
upWAL.OAgene <- upWAL.OA$Row.names

upWAL.RA <- subset(WAL.RA, WAL.RA$Fold.Change >= thresh)
upWAL.RAgene <- upWAL.RA$Row.names


INTUP <- Reduce(intersect, list(upBER.OAgene,upBER.RAgene,upJEN.OAgene,
                                upJEN.RAgene,upBROgene,upVRIgene, upPADgene, 
                                upWAL.RAgene, upWAL.OAgene))



#### DOWN ####
thresh <- -1

downBER.OA <- subset(BER.OA, BER.OA$Fold.Change <= thresh)
downBER.OAgene <- downBER.OA$Gene.Symbol

downBER.RA <- subset(BER.RA, BER.RA$Fold.Change <= thresh)
downBER.RAgene <- downBER.RA$Gene.Symbol

downJEN.OA <- subset(JEN.OA, JEN.OA$Fold.Change <= thresh)
downJEN.OAgene <- downJEN.OA$Gene.Symbol

downJEN.RA <- subset(JEN.RA, JEN.RA$Fold.Change <= thresh)
downJEN.RAgene <- downJEN.RA$Gene.Symbol

downBRO <- subset(BRO, BRO$Fold.Change <= thresh)
downBROgene <- downBRO$Gene.Symbol

downVRI <- subset(VRI, VRI$Fold.Change <= thresh)
downVRIgene <- downVRI$Gene.Symbol

downPAD <- subset(PAD, PAD$Fold.Change <= thresh)
downPADgene <- downPAD$hgnc_symbol

downWAL.OA <- subset(WAL.OA, WAL.OA$Fold.Change <= thresh)
downWAL.OAgene <- downWAL.OA$Row.names

downWAL.RA <- subset(WAL.RA, WAL.RA$Fold.Change <= thresh)
downWAL.RAgene <- downWAL.RA$Row.names

INTDOWN <- Reduce(intersect, list(downBER.OAgene,downBER.RAgene,downJEN.OAgene,
                                  downJEN.RAgene,downBROgene,downVRIgene, downPADgene, 
                                  downWAL.RAgene, downWAL.OAgene))




########################### COMMON GENES ##############################
all <- c(INTUP, INTDOWN)
setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
write.table(INTUP,"upDEGs.txt", row.names = F, col.names = F, quote = F)
write.table(INTDOWN,"downDEGs.txt", row.names = F, col.names = F, quote = F)
write.table(all, "allDEGs.txt", row.names = F, col.names = F, quote = F)




####### Remove Sjogren's #####

setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")

HOR <- read.csv("HOR_filteredresult.csv")
HOR <- HOR[order(HOR$P.Value),]
NAK <- read.csv("NAK_filteredresult.csv")
NAK <- NAK[order(NAK$P.Value),]
MOU <- read.csv("MOU_filteredresult.csv")
MOU <- MOU[order(MOU$P.Value),]



thresh <- 1

upHOR <- subset(HOR, HOR$Fold.Change >= thresh)
upHORgene <- upHOR$Gene.Symbol

upNAK <- subset(NAK, NAK$Fold.Change >= thresh)
upNAKgene <- upNAK$Gene.Symbol

upMOU <- subset(MOU, MOU$Fold.Change >= thresh)
upMOUgene <- upMOU$Gene.Symbol

INTUP_SJO <- Reduce(intersect, list(upHORgene, upNAKgene, upMOUgene))


#### DOWN 
thresh <- -1

downHOR <- subset(HOR, HOR$Fold.Change <= thresh)
downHORgene <- downHOR$Gene.Symbol

downNAK <- subset(NAK, NAK$Fold.Change <= thresh)
downNAKgene <- downNAK$Gene.Symbol

downMOU <- subset(MOU, MOU$Fold.Change <= thresh)
downMOUgene <- downMOU$Gene.Symbol

INTDOWN_SJO <- Reduce(intersect, list(downHORgene, downNAKgene, downMOUgene))


##### COMMON GENES ###
upremove <- Reduce(intersect, list (INTUP, INTUP_SJO))
downremove <- Reduce(intersect, list(INTDOWN, INTDOWN_SJO))


####### Remove SLE #####
setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")

MID <- read.csv("MID_filteredresult.csv")
MID <- MID[order(MID$P.Value),]
SCH <- read.csv("SCH_filteredresult.csv")
SCH <- SCH[order(SCH$P.Value),]

thresh <- 1

upMID <- subset(MID, MID$Fold.Change >= thresh)
upMIDgene <- upMID$Gene.Symbol

upSCH <- subset(SCH, SCH$Fold.Change >= thresh)
upSCHgene <- upSCH$Gene.Symbol

INTUP_SLE <- Reduce(intersect, list(upMIDgene, upSCHgene))


#### DOWN 
thresh <- -1

downMID <- subset(MID, MID$Fold.Change <= thresh)
downMIDgene <- downMID$Gene.Symbol

downSCH <- subset(SCH, SCH$Fold.Change <= thresh)
downSCHgene <- downSCH$Gene.Symbol

INTDOWN_SLE <- Reduce(intersect, list(downMIDgene, downSCHgene))


##### COMMON GENES ###
upremove2 <- Reduce(intersect, list (INTUP, INTUP_SLE))
downremove2 <- Reduce(intersect, list(INTDOWN, INTDOWN_SLE))

# ####### Remove MS #####
# setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
# 
# KEM <- read.csv("KEM_filteredresult.csv")
# KEM <- KEM[order(KEM$P.Value),]
# MAT <- read.csv("MAT_filteredresult.csv")
# MAT <- MAT[order(MAT$P.Value),]
# 
# thresh <- 1
# 
# upKEM <- subset(KEM, KEM$Fold.Change >= thresh)
# upKEMgene <- upKEM$Gene.Symbol
# 
# upMAT <- subset(MAT, MAT$Fold.Change >= thresh)
# upMATgene <- upMAT$Gene.Symbol
# 
# INTUP_MS <- Reduce(intersect, list(upKEMgene, upMATgene))
# 
# 
# #### DOWN 
# thresh <- -1
# 
# downKEM <- subset(KEM, KEM$Fold.Change <= thresh)
# downKEMgene <- downKEM$Gene.Symbol
# 
# downMAT <- subset(MAT, MAT$Fold.Change <= thresh)
# downMATgene <- downMAT$Gene.Symbol
# 
# INTDOWN_MS <- Reduce(intersect, list(downKEMgene, downMATgene))
# 
# 
# ##### COMMON GENES ###
# upremove3 <- Reduce(intersect, list (INTUP, INTUP_MS))
# downremove3 <- Reduce(intersect, list(INTDOWN, INTDOWN_MS))



####### Remove Alopecia #####
setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")

JAB <- read.csv("JAB_filteredresult.csv")
JAB <- JAB[order(JAB$P.Value),]
JAB2 <- read.csv("JAB2_filteredresult.csv")
JAB2 <- JAB2[order(JAB2$P.Value),]

thresh <- 1

upJAB <- subset(JAB, JAB$Fold.Change >= thresh)
upJABgene <- upJAB$Gene.Symbol

upJAB2 <- subset(JAB2, JAB2$Fold.Change >= thresh)
upJAB2gene <- upJAB2$Gene.Symbol

INTUP_AL <- Reduce(intersect, list(upJABgene, upJAB2gene))


#### DOWN
thresh <- -1

downJAB <- subset(JAB, JAB$Fold.Change <= thresh)
downJABgene <- downJAB$Gene.Symbol

downJAB2 <- subset(JAB2, JAB2$Fold.Change <= thresh)
downJAB2gene <- downJAB2$Gene.Symbol

INTDOWN_AL <- Reduce(intersect, list(downJABgene, downJAB2gene))


##### COMMON GENES ###
upremove3 <- Reduce(intersect, list (INTUP, INTUP_AL))
downremove3 <- Reduce(intersect, list(INTDOWN, INTDOWN_AL))


###### REMOVE COMMON GENES ####
# Sjogrens
resultsup <- subset(INTUP, !(INTUP %in% upremove))
resultsdown <- subset(INTDOWN, !(INTDOWN %in% downremove))
results <- c(resultsup, resultsdown)

# SLE
resultsup <- subset(resultsup, !(resultsup %in% upremove2))
resultsdown <- subset(resultsdown, !(resultsdown %in% downremove2))
results <- c(resultsup, resultsdown)

# # MS #ADDITION OF MS FILTER REDUCED SIGNAL FOR RHEUMATOID ARTHRITIS
# resultsup <- subset(resultsup, !(resultsup %in% upremove3))
# resultsdown <- subset(resultsdown, !(resultsdown %in% downremove3))
# results <- c(resultsup, resultsdown)

# Alopecia
resultsup <- subset(resultsup, !(resultsup %in% upremove3))
resultsdown <- subset(resultsdown, !(resultsdown %in% downremove3))
results <- c(resultsup, resultsdown)

setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
write.table(resultsup,"upDEGs_SJO_SLE_ALO.txt", row.names = F, col.names = F, quote = F)
write.table(resultsdown,"downDEGs_SJO_SLE_ALO.txt", row.names = F, col.names = F, quote = F)
write.table(results, "allDEGs_SJO_SLE_ALO.txt", row.names = F, col.names = F, quote = F)


RAGenes <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/RAGenes.txt")
OAGenes <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/OAGenes.txt")

intersect(resultsup, RAGenes)
intersect(resultsup, OAGenes)

intersect(resultsdown, RAGenes)
intersect(resultsdown, OAGenes)

intersect(results, RAGenes)
intersect(results, OAGenes)
