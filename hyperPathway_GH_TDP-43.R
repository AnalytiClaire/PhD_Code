
library(pathprint)
library(hgu133plus2.db)

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GSEA/PCxN Example/probesets/")

### Gene background ###

sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)])
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]
sym.genes <- t(sym.genes)
allgenes <- sym.genes[!duplicated(sym.genes),]

### DEG Thresholds ###

common_DEGs <- read.table(file = "/Users/clairegreen/Desktop/TDP-43DEGs.txt")
common_DEGs <- common_DEGs$V1

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
Y <- read.csv(file = "AllgenesNO.csv", na.strings = c("", "NA)"))
Y <- as.list(Y)
Y <- lapply(Y, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/")
W <- read.csv(file = "BenchmarkGenes.csv", na.strings = c("", "NA)"))
W <- as.list(W)
W<- lapply(W, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint 25.04.16/hyperPathway/All.Pathways/PathprintPathways(29)/")
Z <- read.csv(file = "pathprintgenes.csv", na.strings = c("", "NA)"))
Z <- as.list(Z)
Z <- lapply(Z, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/")
TRL_all <- read.csv(file = "TRL_Comparison.csv", na.strings = c("", "NA)"))
TRL_all <- as.list(Z)
TRL_all <- lapply(Z, function(x) x[!is.na(x)])



### GH Genes ###

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/")

X <- read.table(file = "GH_TRL_1000.txt", header = FALSE)
X <- X$V1


### ALL DEG RESULTS ###

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/SA_Genes/")

SA_TRL <- read.table("tdp43_top_TRL_sigDEG.genes.txt")
SA_CYT <- read.table("tdp43_top_CYT_sigDEG.genes.txt")
SA_WCT <- read.table("tdp43_top_WCT_sigDEG.genes.txt")

SA_TRL_gene <- as.vector(SA_TRL$GeneSymbol)
SA_TRL_gene <- SA_TRL_gene[!is.na(SA_TRL_gene)]
SA_CYT_gene <- as.vector(SA_CYT$GeneSymbol)
SA_CYT_gene <- SA_CYT_gene[!is.na(SA_CYT_gene)]
SA_WCT_gene <- as.vector(SA_WCT$GeneSymbol)
SA_WCT_gene <- SA_WCT_gene[!is.na(SA_TRL_gene)]

### LIMMA ####
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/Limma/")
CG_TRL_limma <- read.csv("GH_HEK_48hHGNC_TRL_rankeduniqueresult.csv")
CG_TRL_sig_limma <- subset(CG_TRL_limma, subset=(P.Value < 0.05))
CG_CYT_limma <- read.csv("GH_HEK_48h_CYT_HGNCrankeduniqueresult.csv")
CG_CYT_sig_limma <- subset(CG_CYT_limma, subset=(adj.P.Val < 0.05))
CG_WCT_limma <- read.csv("GH_HEK_48h_WCT_HGNCrankeduniqueresult.csv")
CG_WCT_sig_limma <- subset(CG_WCT_limma, subset=(adj.P.Val < 0.05))

CG_TRL_gene_limma <- CG_TRL_sig_limma$hgnc_symbol
CG_CYT_gene_limma <- CG_CYT_sig_limma$hgnc_symbol
CG_WCT_gene_limma <- CG_WCT_sig_limma$hgnc_symbol

### DESEQ2 ####
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/DEseq2/")
CG_TRL_deseq2 <- read.csv("TRL_diffexpr-results.csv")
CG_TRL_sig_deseq2 <- subset(CG_TRL_deseq2, subset=(padj < 0.05))
CG_CYT_deseq2 <- read.csv("CYT_diffexpr-results.csv")
CG_CYT_sig_deseq2 <- subset(CG_CYT_deseq2, subset=(padj < 0.05))
CG_WCT_deseq2 <- read.csv("WCT_diffexpr-results.csv")
CG_WCT_sig_deseq2 <- subset(CG_WCT_deseq2, subset=(padj < 0.05))

CG_TRL_gene_deseq2 <- as.vector(CG_TRL_sig_deseq2$hgnc_symbol)
CG_CYT_gene_deseq2 <- as.vector(CG_CYT_sig_deseq2$hgnc_symbol)
CG_WCT_gene_deseq2 <- as.vector(CG_WCT_sig_deseq2$hgnc_symbol)



setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/LC/")
LC_TRL <- read.table("Q331K GRASPS_BitSeq_EdgeR_DE transcripts.txt")
LC_TRL_sig <- subset(LC_TRL, subset=(FDR < 0.05))
LC_CYT <- read.table("Q331K cytoplasmic transcriptomes_BitSeq_EdgeR_DE transcripts.txt")
LC_CYT_sig <- subset(LC_CYT, subset=(FDR < 0.05))
LC_WCT <- read.table("Q331K whole cell transcriptomes_BitSeq_EdgeR_DE transcripts.txt")
LC_WCT_sig <- subset(LC_WCT, subset=(FDR < 0.05))

LC_TRL_gene <- as.vector(LC_TRL_sig$V3)
LC_CYT_gene <- as.vector(LC_CYT_sig$V3)
LC_WCT_gene <- as.vector(LC_WCT_sig$V3)






### Enrichment analysis ###
pathwayEnrichment <- hyperPathway(
  genelist = SA_CYT_gene,
  geneset = W,
  Nchip = length(allgenes)
)
setwd (dir = "/Users/clairegreen/Desktop/")
write.csv(pathwayEnrichment, file = "CG_enrich.csv")




overlap <- Reduce(intersect, list(LC_TRL_gene, CG_TRL_gene_deseq2))
print(overlap)

















# write.table(x = CG_TRL_gene, file = "CG_TRL_gene.txt", quote = FALSE, row.names = FALSE)
# # write.table(x = SA_TRL_gene, file = "SA_TRL_gene.txt", quote = FALSE, row.names = FALSE)
# # write.table(x = LC_TRL_gene, file = "LC_TRL_gene.txt", quote = FALSE, row.names = FALSE)
# # 
# # 
# library(VennDiagram)
# results <- read.csv("TRL_Comparison.csv", na.strings = c("", "NA)"))
# results <- as.list(results)
# # results <- lapply(results, function(x) x[!duplicated(x)])
# results <- lapply(results, function(x) x[!is.na(x)])
# 
# venn <- calculate.overlap(
#   x = list(
#     "SA" = SA_TRL_gene,
#     "LC" = LC_TRL_gene,
#     "CG" = CG_TRL_gene))
# 
# venn <- calculate.overlap(x = results)
# 
# grid.newpage()
# draw.triple.venn(area1 = 697, area2 = 2528, area3 = 2390, n12 = , n23 = , n13 = 35, 
#                  n123 = 0, category = c("Manifesting vs Control", "Non-manifesting vs Control", "Manifesting vs Nonmanifesting"), lty = "blank", 
#                  fill = c("skyblue", "violet", "coral"), cat.dist = -0.1)
# 
# vennDiagram(results)
# 
# ### Intersect ###
# overlap <- Reduce(intersect, list(SA_TRL_gene, W$GeneCards.ALS))
# print(overlap)
# 
# 
# # A simple single-set diagram
# cardiome <- letters[1:10]
# superset <- letters[8:24]
# overlap <- calculate.overlap(
#   x = list(
#     "Cardiome" = cardiome,
#     "SuperSet" = superset
#   )
#);
