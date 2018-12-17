setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
sals <- read.csv("sals_unique.csv")
sals <- sals[order(sals$P.Value),]
ftld <- read.csv("ftld_unique.csv")
ftld <- ftld[order(ftld$P.Value),]
vcp <- read.csv("vcp_unique.csv")
vcp <- vcp[order(vcp$P.Value),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
rav <- read.csv("RAV_results_keepfiltering.csv")


thresh <- 1

upC9 <- subset(C9, C9$Fold.Change >= thresh)
upC9gene <- upC9$Gene.Symbol

upSALS <- subset(sals, sals$Fold.Change >= thresh)
upSALSgene <- upSALS$Gene.Symbol

upFTLD <- subset(ftld, ftld$Fold.Change >= thresh)
upFTLDgene <- upFTLD$Gene.Symbol

upVCP <- subset(vcp, vcp$Fold.Change >= thresh)
upVCPgene <- upVCP$Gene.Symbol

upPET <- subset(pet, pet$FoldChange >= thresh)
upPETgene <- upPET$hgnc_symbol

upRAV <- subset(rav, rav$FoldChange >= thresh)
upRAVgene <- upRAV$hgnc_symbol

INTUP_TDP <- Reduce(intersect, list(upC9gene, upSALSgene, upFTLDgene, upVCPgene, upPETgene, upRAVgene))

# cat(INTUP, sep = "\n")

# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults")
# write.table(INTUP, "intersect_up_1.txt", col.names = F, row.names = F, quote = F)



#### DOWN ####
thresh <- -1

downC9 <- subset(C9, C9$Fold.Change <= thresh)
downC9gene <- downC9$Gene.Symbol

downSALS <- subset(sals, sals$Fold.Change <= thresh)
downSALSgene <- downSALS$Gene.Symbol

downFTLD <- subset(ftld, ftld$Fold.Change <= thresh)
downFTLDgene <- downFTLD$Gene.Symbol

downVCP <- subset(vcp, vcp$Fold.Change <= thresh)
downVCPgene <- downVCP$Gene.Symbol

downPET <- subset(pet, pet$FoldChange <= thresh)
downPETgene <- downPET$hgnc_symbol

downRAV <- subset(rav, rav$FoldChange <= thresh)
downRAVgene <- downRAV$hgnc_symbol

INTDOWN_TDP <- Reduce(intersect, list(downC9gene, downSALSgene, downFTLDgene, downVCPgene, downPETgene, downRAVgene))

# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults")
# write.table(INTDOWN, "intersect_down_1.txt", col.names = F, row.names = F, quote = F)

# cat(INTDOWN, sep = "\n")



############################################################################
################### NON-TDP ALS ############################################
############################################################################

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/non-TDP/")

FUS <- read.csv("FUSrankeduniqueresult.csv")
FUS <- FUS[order(FUS$P.Value),]
SOD1 <- read.csv("SOD1rankeduniqueresult.csv")
SOD1 <- SOD1[order(SOD1$P.Value),]
# CBFTLD <- read.csv("CBFTLDrankeduniqueresult.csv")
# CBFTLD <- CBFTLD[order(CBFTLD$P.Value),]
# FIB <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/PH_Fibroblasts/WCT/WCT_norm.csv")
# FIB <- FIB[order(FIB$pvalue),]

#### UP ####
thresh <- 1

upFUS <- subset(FUS, FUS$Fold.Change >= thresh)
upFUSgene <- upFUS$Gene.Symbol

upSOD1 <- subset(SOD1, SOD1$Fold.Change >= thresh)
upSOD1gene <- upSOD1$Gene.Symbol

INTUP_non <- Reduce(intersect, list(upSOD1gene, upFUSgene))

#### DOWN ####
thresh <- -1

downFUS <- subset(FUS, FUS$Fold.Change <= thresh)
downFUSgene <- downFUS$Gene.Symbol

downSOD1 <- subset(SOD1, SOD1$Fold.Change <= thresh)
downSOD1gene <- downSOD1$Gene.Symbol

INTDOWN_non <- Reduce(intersect, list(downSOD1gene, downFUSgene))



########################### COMMON GENES ##############################
upremove <- Reduce(intersect, list (INTUP_TDP, INTUP_non))
downremove <- Reduce(intersect, list(INTDOWN_TDP, INTDOWN_non))



################ REMOVE COMMON GENES ######################
resultsup <- subset(INTUP_TDP, !(INTUP_TDP %in% upremove))
resultsdown <- subset(INTDOWN_TDP, !(INTDOWN_TDP %in% downremove))
results <- c(resultsup, resultsdown)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/CurrentStuff")
write.table(resultsup, "Filtered_up_genes.txt", quote = F, row.names = F, col.names = F)
write.table(resultsdown, "Filtered_down_genes.txt", quote = F, row.names = F, col.names = F)
write.table(results, "Filtered_all_genes.txt", quote = F, row.names = F, col.names = F)
# cat(resultsup, sep="\n")
