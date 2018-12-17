setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

thresh <- 5996

C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
C96500 <- C9[1:thresh,]
CH <- read.csv("CH_unique.csv")
CH <- CH[order(CH$P.Value),]
CH6500 <- CH[1:thresh,]
sals <- read.csv("sals_unique.csv")
sals <- sals[order(sals$P.Value),]
sals6500 <- sals[1:thresh,]
ftld <- read.csv("ftld_unique.csv")
ftld <- ftld[order(ftld$P.Value),]
ftld6500 <- ftld[1:thresh,]
vcp <- read.csv("vcp_unique.csv")
vcp <- vcp[order(vcp$P.Value),]
vcp6500 <- vcp[1:thresh,]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_results_keepfiltering.csv")
pet <- pet[!duplicated(pet$hgnc_symbol),]
pet6500 <- pet[1:thresh,]
rav <- read.csv("RAV_results_keepfiltering.csv")
rav <- rav[!duplicated(rav$hgnc_symbol),]
rav6500 <- rav[1:thresh,]



C9up <-subset(C96500, subset=(logFC > 0))
C9upgene <- C9up$Gene.Symbol
C9down <-subset(C96500, subset=(logFC < 0))
C9downgene <- C9down$Gene.Symbol

CHup <-subset(CH6500, subset=(logFC > 0))
CHupgene <- CHup$Gene.Symbol
CHdown <-subset(CH6500, subset=(logFC < 0))
CHdowngene <- CHdown$Gene.Symbol

salsup <-subset(sals6500, subset=(logFC > 0))
salsupgene <- salsup$Gene.Symbol
salsdown <-subset(sals6500, subset=(logFC < 0))
salsdowngene <- salsdown$Gene.Symbol

ftldup <-subset(ftld6500, subset=(logFC > 0))
ftldupgene <- ftldup$Gene.Symbol
ftlddown <-subset(ftld6500, subset=(logFC < 0))
ftlddowngene <- ftlddown$Gene.Symbol

vcpup <-subset(vcp6500, subset=(logFC > 0))
vcpupgene <- vcpup$Gene.Symbol
vcpdown <-subset(vcp6500, subset=(logFC < 0))
vcpdowngene <- vcpdown$Gene.Symbol

petup <-subset(pet6500, subset=(log2FoldChange > 0))
petupgene <- petup$hgnc_symbol
petdown <-subset(pet6500, subset=(log2FoldChange < 0))
petdowngene <- petdown$hgnc_symbol

ravup <-subset(rav6500, subset=(log2FoldChange > 0))
ravupgene <- ravup$hgnc_symbol
ravdown <-subset(rav6500, subset=(log2FoldChange < 0))
ravdowngene <- ravdown$hgnc_symbol

intersect_up <- Reduce(intersect, list(C9upgene, CHupgene, salsupgene, ftldupgene, vcpupgene, petupgene, ravupgene ))
cat(intersect_up, sep = "\n")
intersect_down <- Reduce(intersect, list(C9downgene, CHdowngene, salsdowngene, ftlddowngene, vcpdowngene, petdowngene, ravdowngene ))
cat(intersect_down, sep = "\n")

write.table(intersect_up, "intersect_up.txt", quote = F, col.names = F, row.names = F)
write.table(intersect_down, "intersect_down.txt", quote = F, col.names = F, row.names = F)
