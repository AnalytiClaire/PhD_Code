#############################################################
################### FOLD CHANGE DEG #########################
#############################################################

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

LEW <- read.csv("LEWfilteredresult.csv")
MID3 <- read.csv("MID3filteredresult.csv")
MID4 <- read.csv("MID4filteredresult.csv")
MOR.FC <- read.csv("MOR.FCfilteredresult.csv")
DIJ <- read.csv("DIJfilteredresult.csv")
FFR <- read.csv("FFRfilteredresult.csv")
MID1 <- read.csv("MID1filteredresult.csv")
MID2 <- read.csv("MID2filteredresult.csv")
MOR.SN <- read.csv("MOR.SNfilteredresult.csv")
DUM <- read.csv("DUM_UniqueGene_DESeq2.csv")
BOT <- read.csv("BOTrankeduniqueresult.csv")
BOT2 <- read.csv("BOT2rankeduniqueresult.csv")




thresh <- 1

upLEW <- subset(LEW, LEW$Fold.Change >= thresh)
upLEWgene <- upLEW$Gene.Symbol

upMID3<- subset(MID3, MID3$Fold.Change >= thresh)
upMID3gene <- upMID3$Gene.Symbol

upMID4 <- subset(MID4, MID4$Fold.Change >= thresh)
upMID4gene <- upMID4$Gene.Symbol

upMOR.FC <- subset(MOR.FC, MOR.FC$Fold.Change >= thresh)
upMOR.FCgene <- upMOR.FC$Gene.Symbol

upDUM <- subset(DUM, DUM$log2FoldChange >= 0)
upDUMgene <- upDUM$hgnc_symbol

upDIJ <- subset(DIJ, DIJ$Fold.Change >= thresh)
upDIJgene <- upDIJ$Gene.Symbol

upFFR<- subset(FFR, FFR$Fold.Change >= thresh)
upFFRgene <- upFFR$Gene.Symbol

upMID1<- subset(MID1, MID1$Fold.Change >= thresh)
upMID1gene <- upMID1$Gene.Symbol

upMID2 <- subset(MID2, MID2$Fold.Change >= thresh)
upMID2gene <- upMID2$Gene.Symbol

upMOR.SN <- subset(MOR.SN, MOR.SN$Fold.Change >= thresh)
upMOR.SNgene <- upMOR.SN$Gene.Symbol

upBOT <- subset(BOT, BOT$Fold.Change >= thresh)
upBOTgene <- upBOT$Gene.Symbol

upBOT2 <- subset(BOT2, BOT2$Fold.Change >= thresh)
upBOT2gene <- upBOT2$Gene.Symbol


INTUP <- Reduce(intersect, list(upLEWgene, upMID3gene, upMID4gene, upMOR.FCgene,
                                upDIJgene, upFFRgene, upMID1gene, upMID2gene, 
                                upMOR.SNgene, upDUMgene, upBOTgene, upBOT2gene))



#### DOWN ####
thresh <- -1

downLEW <- subset(LEW, LEW$Fold.Change <= thresh)
downLEWgene <- downLEW$Gene.Symbol

downMID3 <- subset(MID3, MID3$Fold.Change <= thresh)
downMID3gene <- downMID3$Gene.Symbol

downMID4 <- subset(MID4, MID4$Fold.Change <= thresh)
downMID4gene <- downMID4$Gene.Symbol

downMOR.FC <- subset(MOR.FC, MOR.FC$Fold.Change <= thresh)
downMOR.FCgene <- downMOR.FC$Gene.Symbol

downDUM <- subset(DUM, DUM$log2FoldChange <= 0)
downDUMgene <- downDUM$hgnc_symbol

downDIJ <- subset(DIJ, DIJ$Fold.Change <= thresh)
downDIJgene <- downDIJ$Gene.Symbol

downFFR<- subset(FFR, FFR$Fold.Change <= thresh)
downFFRgene <- downFFR$Gene.Symbol

downMID1<- subset(MID1, MID1$Fold.Change <= thresh)
downMID1gene <- downMID1$Gene.Symbol

downMID2 <- subset(MID2, MID2$Fold.Change <= thresh)
downMID2gene <- downMID2$Gene.Symbol

downMOR.SN <- subset(MOR.SN, MOR.SN$Fold.Change <= thresh)
downMOR.SNgene <- downMOR.SN$Gene.Symbol


downBOT <- subset(BOT, BOT$Fold.Change <= thresh)
downBOTgene <- downBOT$Gene.Symbol

downBOT2 <- subset(BOT2, BOT2$Fold.Change <= thresh)
downBOT2gene <- downBOT2$Gene.Symbol

INTDOWN <- Reduce(intersect, list(downLEWgene, downMID3gene, downMID4gene, downMOR.FCgene,
                                  downDIJgene, downFFRgene, downMID1gene, downMID2gene, 
                                  downMOR.SNgene, downDUMgene, downBOTgene, downBOT2gene))






####### ALS signature ##########

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
sals <- read.csv("sals_unique.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
rav <- read.csv("RAV_results_keepfiltering.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/non-TDP/")

FUS <- read.csv("FUSrankeduniqueresult.csv")
SOD1 <- read.csv("SOD1rankeduniqueresult.csv")


thresh <- 1

upC9 <- subset(C9, C9$Fold.Change >= thresh)
upC9gene <- upC9$Gene.Symbol

upSALS <- subset(sals, sals$Fold.Change >= thresh)
upSALSgene <- upSALS$Gene.Symbol

upPET <- subset(pet, pet$FoldChange >= thresh)
upPETgene <- upPET$hgnc_symbol

upRAV <- subset(rav, rav$FoldChange >= thresh)
upRAVgene <- upRAV$hgnc_symbol

upFUS <- subset(FUS, FUS$Fold.Change >= thresh)
upFUSgene <- upFUS$Gene.Symbol

upSOD1 <- subset(SOD1, SOD1$Fold.Change >= thresh)
upSOD1gene <- upSOD1$Gene.Symbol

INTUP_ALS <- Reduce(intersect, list(upC9gene, upSALSgene, upPETgene, upRAVgene, upFUSgene, upSOD1gene))


#### DOWN 
thresh <- -1

downC9 <- subset(C9, C9$Fold.Change <= thresh)
downC9gene <- downC9$Gene.Symbol

downSALS <- subset(sals, sals$Fold.Change <= thresh)
downSALSgene <- downSALS$Gene.Symbol

downPET <- subset(pet, pet$FoldChange <= thresh)
downPETgene <- downPET$hgnc_symbol

downRAV <- subset(rav, rav$FoldChange <= thresh)
downRAVgene <- downRAV$hgnc_symbol

downFUS <- subset(FUS, FUS$Fold.Change <= thresh)
downFUSgene <- downFUS$Gene.Symbol

downSOD1 <- subset(SOD1, SOD1$Fold.Change <= thresh)
downSOD1gene <- downSOD1$Gene.Symbol

INTDOWN_ALS <- Reduce(intersect, list(downC9gene, downSALSgene, downPETgene, downRAVgene, downFUSgene, downSOD1gene))





####### sPD Blood signature ##########

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

AMA <- read.csv("AMAfilteredresult.csv")
RON <- read.csv("RONfilteredresult.csv")


thresh <- 1

upAMA <- subset(AMA, AMA$Fold.Change >= thresh)
upAMAgene <- upAMA$Gene.Symbol

upRON <- subset(RON, RON$Fold.Change >= thresh)
upRONgene <- upRON$Gene.Symbol


INTUP_blood <- Reduce(intersect, list(upAMAgene, upRONgene))


#### DOWN ###
thresh <- -1

downAMA <- subset(AMA, AMA$Fold.Change <= thresh)
downAMAgene <- downAMA$Gene.Symbol

downRON <- subset(RON, RON$Fold.Change <= thresh)
downRONgene <- downRON$Gene.Symbol


INTDOWN_blood <- Reduce(intersect, list(downAMAgene, downRONgene))






####### fam Blood signature ##########

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

AMA_PINK1 <- read.csv("AMA_PINK1filteredresult.csv")
AMA_PARKIN <- read.csv("AMA_PARKINfilteredresult.csv")
AMA_ATP13A2 <- read.csv("AMA_ATP13A2filteredresult.csv")


thresh <- 1

upAMA_PINK1 <- subset(AMA_PINK1, AMA_PINK1$Fold.Change >= thresh)
upAMA_PINK1gene <- upAMA_PINK1$Gene.Symbol

upAMA_PARKIN <- subset(AMA_PARKIN, AMA_PARKIN$Fold.Change >= thresh)
upAMA_PARKINgene <- upAMA_PARKIN$Gene.Symbol

upAMA_ATP13A2 <- subset(AMA_ATP13A2, AMA_ATP13A2$Fold.Change >= thresh)
upAMA_ATP13A2gene <- upAMA_ATP13A2$Gene.Symbol


INTUP_famblood <- Reduce(intersect, list(upAMA_PINK1gene, upAMA_PARKINgene, upAMA_ATP13A2gene))


#### DOWN ###
thresh <- -1

downAMA_PINK1 <- subset(AMA_PINK1, AMA_PINK1$Fold.Change <= thresh)
downAMA_PINK1gene <- downAMA_PINK1$Gene.Symbol

downAMA_PARKIN <- subset(AMA_PARKIN, AMA_PARKIN$Fold.Change <= thresh)
downAMA_PARKINgene <- downAMA_PARKIN$Gene.Symbol

downAMA_ATP13A2 <- subset(AMA_ATP13A2, AMA_ATP13A2$Fold.Change <= thresh)
downAMA_ATP13A2gene <- downAMA_ATP13A2$Gene.Symbol


INTDOWN_famblood <- Reduce(intersect, list(downAMA_PINK1gene, downAMA_PARKINgene, downAMA_ATP13A2gene))





################
### REMOVALS ###
################

##### COMMON GENES ###
upremove <- Reduce(intersect, list (INTUP, INTUP_ALS))
downremove <- Reduce(intersect, list(INTDOWN, INTDOWN_ALS))

##### REMOVE COMMON GENES ###
resultsup <- subset(INTUP, !(INTUP %in% upremove))
resultsdown <- subset(INTDOWN, !(INTDOWN %in% downremove))
results <- c(resultsup, resultsdown)

##### COMMON GENES ###
upremove2 <- Reduce(intersect, list(resultsup, INTUP_blood))
downremove2 <- Reduce(intersect, list(resultsdown, INTDOWN_blood))

##### REMOVE COMMON GENES ###
resultsup <- subset(resultsup, !(resultsup %in% upremove2))
resultsdown <- subset(resultsdown, !(resultsdown %in% downremove2))
results <- c(resultsup, resultsdown)

##### COMMON GENES ###
upremove3 <- Reduce(intersect, list(resultsup, INTUP_famblood))
downremove3 <- Reduce(intersect, list(resultsdown, INTDOWN_famblood))

##### REMOVE COMMON GENES ###
resultsup <- subset(resultsup, !(resultsup %in% upremove3))
resultsdown <- subset(resultsdown, !(resultsdown %in% downremove3))
results <- c(resultsup, resultsdown)






setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/")
write.table(resultsup, "ALS_sfblood_UPgenes.txt", quote = F, row.names = F, col.names = F)
write.table(resultsdown, "ALS_sfblood_DOWNgenes.txt", quote = F, row.names = F, col.names = F)
write.table(results, "ALS_sfblood_ALLgenes.txt", quote = F, row.names = F, col.names = F)



PDGenes <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/ParkinsonsDiseaseMalacards.txt")
LRRK2PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/ReyniersLRRK2PPI.txt")
LRRK2BioGrid <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/BioGrid_LRRK2.txt")
DEGs <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ALS_sfblood_ALLgenes.txt")

intersect(DEGs, PDGenes)
intersect(DEGs, LRRK2PPI)
intersect(DEGs, LRRK2BioGrid)

