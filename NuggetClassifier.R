setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression")

C9 <- read.csv("C9rankeduniqueresult.csv")
sals <- read.csv("sALSrankeduniqueresult.csv")
ftld <- read.csv("FTLDrankeduniqueresult.csv")
vcp <- read.csv("VCPrankeduniqueresult.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
pet <- pet[!duplicated(pet$hgnc_symbol),]
rav <- read.csv("RAV_results_keepfiltering.csv")
rav <- rav[!duplicated(rav$hgnc_symbol),]



#Find genes that are present across all FULL datasets
C9gen <- C9$Gene.Symbol
salsgen <- sals$Gene.Symbol
ftldgen <- ftld$Gene.Symbol
vcpgen <- vcp$Gene.Symbol
petgen <- pet$hgnc_symbol
ravgen <- rav$hgnc_symbol


genelist <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NuggetGenes.txt")


subsetC9 <- subset(C9, C9$Gene.Symbol %in% genelist, drop = TRUE)
subsetC9 <- subsetC9[!duplicated(subsetC9$Gene.Symbol),]
rownames(subsetC9) <- subsetC9$Gene.Symbol
subsetC9[,1] = NULL
subsetC9 <- subsetC9[order(row.names(subsetC9)),]

subsetsals <- subset(sals, sals$Gene.Symbol %in% genelist, drop = TRUE)
subsetsals <- subsetsals[!duplicated(subsetsals$Gene.Symbol),]
rownames(subsetsals) <- subsetsals$Gene.Symbol
subsetsals[,1] = NULL
subsetsals <- subsetsals[order(row.names(subsetsals)),]

subsetftld <- subset(ftld, ftld$Gene.Symbol %in% genelist, drop = TRUE)
subsetftld <- subsetftld[!duplicated(subsetftld$Gene.Symbol),]
rownames(subsetftld) <- subsetftld$Gene.Symbol
subsetftld[,1] = NULL
subsetftld <- subsetftld[order(row.names(subsetftld)),]

subsetvcp <- subset(vcp, vcp$Gene.Symbol %in% genelist, drop = TRUE)
subsetvcp <- subsetvcp[!duplicated(subsetvcp$Gene.Symbol),]
rownames(subsetvcp) <- subsetvcp$Gene.Symbol
subsetvcp[,1] = NULL
subsetvcp <- subsetvcp[order(row.names(subsetvcp)),]

subsetpet <- subset(pet, pet$hgnc_symbol %in% genelist, drop = TRUE)
subsetpet <- subsetpet[!duplicated(subsetpet$hgnc_symbol),]
rownames(subsetpet) <- subsetpet$Gene.Symbol
subsetpet[,1] = NULL
subsetpet <- subsetpet[order(subsetpet$hgnc_symbol),]
rownames(subsetpet) <- subsetpet$hgnc_symbol

subsetrav <- subset(rav, rav$hgnc_symbol %in% genelist, drop = TRUE)
subsetrav <- subsetrav[!duplicated(subsetrav$hgnc_symbol),]
rownames(subsetrav) <- subsetrav$Gene.Symbol
subsetrav[,1] = NULL
subsetrav <- subsetrav[order(subsetrav$hgnc_symbol),]
rownames(subsetrav) <- subsetrav$hgnc_symbol



intersect <- Reduce(intersect, list(rownames(subsetC9), rownames(subsetsals), rownames(subsetftld), rownames(subsetvcp), 
                                   subsetpet$hgnc_symbol, subsetrav$hgnc_symbol))


#Generate a matrix with gene names and log fold change values from each
LFC <- data.frame(gene=row.names(subsetC9), 
                  C9orf72=subsetC9$logFC,
                  sALS.1=subsetsals$logFC, 
                  FTLD=subsetftld$logFC,
                  VCP=subsetvcp$logFC,
                  PET=subsetpet$log2FoldChange, 
                  RAV=subsetrav$log2FoldChange)

rownames(LFC) <- LFC$gene
LFC[,1] <- NULL
LFC <- as.matrix(LFC)
hmcols<- rev(brewer.pal(9, "RdYlBu"))
heatmap(LFC)


Patvcon.C9 <- subsetC9[,20:30]
colnames(Patvcon.C9) <- c(rep("Control", 3), rep("Patient", 8))
Patvcon.C9 <- as.matrix(Patvcon.C9)
heatmap(Patvcon.C9, 
        cexCol = 1, 
        main = "C9orf72")

Patvcon.sals <- subsetsals[,20:29]
colnames(Patvcon.sals) <- c(rep("Control", 3), rep("Patient", 7))
Patvcon.sals <- as.matrix(Patvcon.sals)
heatmap(Patvcon.sals, 
        cexCol = 1, 
        main = "sALS")

Patvcon.ftld <- subsetftld[,19:42]
colnames(Patvcon.ftld) <- c(rep("Control", 8), rep("Patient", 16))
Patvcon.ftld <- as.matrix(Patvcon.ftld)
heatmap(Patvcon.ftld, 
        cexCol = 1, 
        main = "FTLD")

Patvcon.vcp <- subsetvcp[,20:29]
colnames(Patvcon.vcp) <- c(rep("Control", 3), rep("Patient", 7))
Patvcon.vcp <- as.matrix(Patvcon.vcp)
heatmap(Patvcon.vcp, 
        cexCol = 1, 
        main = "VCP")

Patvcon.pet <- subsetpet[,9:34]
colnames(Patvcon.pet) <- c(rep("Control", 9), rep("Patient", 17))
Patvcon.pet <- as.matrix(Patvcon.pet)
heatmap(Patvcon.pet, 
        cexCol = 1, 
        main = "PET")

Patvcon.rav <- subsetrav[,9:29]
colnames(Patvcon.rav) <- c(rep("Control", 8), rep("Patient", 13))
Patvcon.rav <- as.matrix(Patvcon.rav)
heatmap(Patvcon.rav,
        cexCol = 1, 
        main = "RAV")




setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/non-TDP/ALL/")
fus <- read.csv("FUSrankeduniqueresult.csv")
sod1 <- read.csv("SOD1rankeduniqueresult.csv")
# cbftld <- read.csv("CBFTLDrankeduniqueresult.csv")
phfib <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/PH_Fibroblasts/WCT_Alltranscripts_norm.csv")

subsetphfib <- subset(phfib, phfib$hgnc_symbol %in% genelist, drop = TRUE)
subsetphfib <- subsetphfib[!duplicated(subsetphfib$hgnc_symbol),]
rownames(subsetphfib) <- subsetphfib$hgnc_symbol
subsetphfib <- subsetphfib[order(subsetphfib$hgnc_symbol),]

subsetfus <- subset(fus, fus$Gene.Symbol %in% genelist, drop = TRUE)
subsetfus <- subsetfus[!duplicated(subsetfus$Gene.Symbol),]
rownames(subsetfus) <- subsetfus$Gene.Symbol
subsetfus <- subsetfus[order(subsetfus$Gene.Symbol),]

subsetsod1 <- subset(sod1, sod1$Gene.Symbol %in% genelist, drop = TRUE)
subsetsod1 <- subsetsod1[!duplicated(subsetsod1$Gene.Symbol),]
rownames(subsetsod1) <- subsetsod1$Gene.Symbol
subsetsod1 <- subsetsod1[order(subsetsod1$Gene.Symbol),]





LFC <- data.frame(row.names =row.names(subsetC9), 
                  C9orf72=subsetC9$logFC,
                  sALS=subsetsals$logFC, 
                  FTLD=subsetftld$logFC,
                  VCP=subsetvcp$logFC, 
                  PET=subsetpet$log2FoldChange, 
                  RAV=subsetrav$log2FoldChange, 
                  PHFib = subsetphfib$log2FoldChange, 
                  FUS = subsetfus$logFC,
                  SOD1 = subsetsod1$logFC)
LFC <- as.matrix(LFC)

FC <- data.frame(row.names =row.names(subsetC9), 
                  C9orf72=subsetC9$Fold.Change,
                  sALS=subsetsals$Fold.Change, 
                  FTLD=subsetftld$Fold.Change,
                  VCP=subsetvcp$Fold.Change, 
                  PET=subsetpet$FoldChange, 
                  RAV=subsetrav$FoldChange, 
                  FUS = subsetfus$Fold.Change,
                  SOD1 = subsetsod1$Fold.Change)



library(gplots)
library(heatmap3)
library(RColorBrewer)
hmcols<- rev(brewer.pal(9, "RdYlBu"))


cols <- c(rep("magenta4",6),
                    rep("magenta2", 3),
                    "magenta4")

heatmap(LFC,
         cexRow = 0.7)





