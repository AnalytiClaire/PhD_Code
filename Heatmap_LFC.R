setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
# CH <- read.csv("CH_unique.csv")
sals <- read.csv("sals_unique.csv")
ftld <- read.csv("ftld_unique.csv")
vcp <- read.csv("vcp_unique.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
pet <- pet[!duplicated(pet$hgnc_symbol),]
rav <- read.csv("RAV_results_keepfiltering.csv")
rav <- rav[!duplicated(rav$hgnc_symbol),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/non-TDP/")
fus <- read.csv("FUSrankeduniqueresult.csv")
sod1 <- read.csv("SOD1rankeduniqueresult.csv")
cbftld <- read.csv("CBFTLDrankeduniqueresult.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/EL_TDP-43/")
EL <- read.csv("EL_results_19052017_normalised_pairwise.csv")



#Find genes that are present across all FULL datasets
C9gen <- C9$Gene.Symbol
# CHgen <- CH$Gene.Symbol
salsgen <- sals$Gene.Symbol
ftldgen <- ftld$Gene.Symbol
vcpgen <- vcp$Gene.Symbol
petgen <- pet$hgnc_symbol
ravgen <- rav$hgnc_symbol
fusgen <- fus$Gene.Symbol
sod1gen <- sod1$Gene.Symbol
cbftldgen <- cbftld$Gene.Symbol
elgen <- EL$hgnc_symbol

genelist <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_upanddown.txt")
genelist <- genelist$V1


#take that subset from each of the datasets

# 
# subsetCH <- subset(CH, CH$Gene.Symbol %in% genelist, drop = TRUE)
# rownames(subsetCH) <- subsetCH$Gene.Symbol
# subsetCH[,1] = NULL
# subsetCH <- subsetCH[order(row.names(subsetCH)),]

subsetC9 <- subset(C9, C9$Gene.Symbol %in% genelist, drop = TRUE)
rownames(subsetC9) <- subsetC9$Gene.Symbol
subsetC9[,1] = NULL
subsetC9 <- subsetC9[order(row.names(subsetC9)),]

subsetsals <- subset(sals, sals$Gene.Symbol %in% genelist, drop = TRUE)
rownames(subsetsals) <- subsetsals$Gene.Symbol
subsetsals[,1] = NULL
subsetsals <- subsetsals[order(row.names(subsetsals)),]

subsetftld <- subset(ftld, ftld$Gene.Symbol %in% genelist, drop = TRUE)
rownames(subsetftld) <- subsetftld$Gene.Symbol
subsetftld[,1] = NULL
subsetftld <- subsetftld[order(row.names(subsetftld)),]

subsetvcp <- subset(vcp, vcp$Gene.Symbol %in% genelist, drop = TRUE)
rownames(subsetvcp) <- subsetvcp$Gene.Symbol
subsetvcp[,1] = NULL
subsetvcp <- subsetvcp[order(row.names(subsetvcp)),]

subsetpet <- subset(pet, pet$hgnc_symbol %in% genelist, drop = TRUE)
rownames(subsetpet) <- subsetpet$Gene.Symbol
subsetpet[,1] = NULL
subsetpet <- subsetpet[order(subsetpet$hgnc_symbol),]

subsetrav <- subset(rav, rav$hgnc_symbol %in% genelist, drop = TRUE)
rownames(subsetrav) <- subsetrav$Gene.Symbol
subsetrav[,1] = NULL
subsetrav <- subsetrav[order(subsetrav$hgnc_symbol),]

subsetsod1 <- subset(sod1, sod1$Gene.Symbol %in% genelist, drop = TRUE)
rownames(subsetsod1) <- subsetsod1$Gene.Symbol
subsetsod1[,1] = NULL
subsetsod1 <- subsetsod1[order(subsetsod1$Gene.Symbol),]

subsetfus <- subset(fus, fus$Gene.Symbol %in% genelist, drop = TRUE)
rownames(subsetfus) <- subsetfus$Gene.Symbol
subsetfus[,1] = NULL
subsetfus <- subsetfus[order(subsetfus$Gene.Symbol),]

subsetcbftld <- subset(cbftld, cbftld$Gene.Symbol %in% genelist, drop = TRUE)
rownames(subsetcbftld) <- subsetcbftld$Gene.Symbol
subsetcbftld[,1] = NULL
subsetcbftld <- subsetcbftld[order(subsetcbftld$Gene.Symbol),]

subsetel<- subset(EL, EL$hgnc_symbol %in% genelist, drop = TRUE)
rownames(subsetel) <- subsetel$hgnc_symbol
subsetel[,1] = NULL
subsetel <- subsetel[order(subsetel$hgnc_symbol),]




genelist <- Reduce(intersect, list(rownames(subsetC9), rownames(subsetsals), rownames(subsetftld), rownames(subsetvcp), 
                                   subsetpet$hgnc_symbol, subsetrav$hgnc_symbol, rownames(subsetsod1), rownames(subsetfus),
                                   rownames(subsetcbftld), rownames(subsetel)))

##### GO ROUND AGAIN #### run genelist through to end and then repeat subsetting

#Generate a matrix with gene names and log fold change values from each
LFC <- data.frame(gene=row.names(subsetC9), 
                C9orf72=subsetC9$logFC,
                sALS.1=subsetsals$logFC, 
                FTLD=subsetftld$logFC,
                VCP=subsetvcp$logFC, 
                C9sALS=subsetpet$log2FoldChange, 
                sALS.2=subsetrav$log2FoldChange,
                FUS=subsetfus$logFC,
                SOD1=subsetsod1$logFC,
                CBFTLD=subsetcbftld$logFC,
                EL=subsetel$log2FoldChange)
                
rownames(LFC) <- LFC$gene
LFC[,1] <- NULL
LFC <- as.matrix(LFC)

# LFCNC <- data.frame(gene=row.names(subsetC9), 
#                    C9orf72=subsetC9$logFC,
#                    sALS=subsetsals$logFC, 
#                    FTLD=subsetftld$logFC,
#                    VCP=subsetvcp$logFC, 
#                    C9_sALS=subsetpet$log2FoldChange, 
#                    sALS=subsetrav$log2FoldChange) 
# 
# rownames(LFCNC) <- LFCNC$gene
# LFCNC[,1] <- NULL
# LFCNC <- as.matrix(LFCNC)

library(gplots)
library(heatmap3)
hmcols<- colorRampPalette(c("green4","green","white", "red","red4"))(256)
# 
# Platform <- c("orange","orange","orange","orange","blue", "blue")
# Tissue <- c("magenta", "magenta", "yellow", "cyan", "yellow", "magenta")

# colsidecolorsLFCNC <- cbind(Platform=c("darkgoldenrod3","darkgoldenrod3","darkgoldenrod3","darkgoldenrod3","darkgoldenrod1", "darkgoldenrod1"),
#                        Disease=c("dodgerblue4","dodgerblue4","cyan","cornflowerblue","dodgerblue4","dodgerblue4"),
#                        Tissue=c("magenta", "magenta", "purple", "purple4", "purple", "magenta"))

colsidecolorsLFC <- cbind(Disease=c(rep("skyblue4",2),
                                    "skyblue3",
                                    "skyblue1",
                                    rep("skyblue4",4),
                                    "skyblue3",
                                    "skyblue4"),
                          Platform=c("goldenrod3","goldenrod3","goldenrod3","goldenrod3",
                                     "gold2", "gold2",
                                     "goldenrod3","goldenrod3","goldenrod3",
                                     "gold2"),
                          Tissue=c(rep("gray25",2),
                                   "gray60",
                                   "gray80",
                                   rep("gray25",4),
                                   "gray95",
                                   "gray80"),
                          Pathology=c(rep("magenta4",6),
                                      rep("magenta2", 3),
                                      "magenta4"))


heatmap3(LFC, 
         col = hmcols,
         Rowv = TRUE,
         Colv = TRUE,
         distfun = dist,
         # hclustfun = hclust,
         scale = "row",
         labCol = colnames(LFC),
         ColSideColors = colsidecolorsLFC,
         ColSideWidth = 1,
         cexCol = 1.3,
         cexRow = 0.01,
         margins = c(7,7),
         xlab = "Datasets",
         ylab = "Genes",
         legendfun = function() showLegend(legend = c("TDP positive", "TDP negative",
                                                      "Microarray", "RNA-seq",
                                                      "Spinal cord", "Frontal Cortex", "Muscle", "Cerebellum",
                                                      "ALS", "FTLD", "IBMPFTD"), 
                                           col = c("magenta4","magenta2",
                                                   "goldenrod3", "gold2",
                                                   "gray25","gray60","gray80","gray95",
                                                   "skyblue4","skyblue3","skyblue1")))
























# heatmap3(LFCNC, 
#          col = hmcols,
#          Rowv = TRUE,
#          Colv = TRUE,
#          distfun = dist,
#          hclustfun = hclust,
#          scale = "row",
#          labCol = colnames(LFCNC),
#          ColSideColors = colsidecolorsLFCNC,
#          ColSideWidth = 1,
#          cexCol = 1.2,
#          cexRow = 0.01,
#          xlab = "Datasets",
#          ylab = "Genes")

# heatmap.2(LFCNC, 
#           Rowv = TRUE, 
#           Colv = TRUE, 
#           scale = "row",
#           distfun = dist, 
#           hclustfun = hclust, 
#           dendrogram = "both",
#           col=hmcols, 
#           ColSideColors = Platform,
#           key = TRUE, 
#           cexCol = 1.2,
#           cexRow = 0.01,
#           srtCol = 45,
#           trace = "none")
# 
# 
# 
# library(heatmap.plus)
# heatmap.plus(LFCNC, 
#              Rowv = TRUE, 
#              Colv = NA, 
#              scale = "row",
#              distfun = dist, 
#              col = hmcols,
#              hclustfun = hclust,
#              cexCol = 1.2,
#              cexRow = 0.01,
#              ColSideColors = colsidecolors)


# 
# Disease = HeatmapAnnotation(df = data.frame(Condition = c(rep("C9", 1), 
#                                                           rep("CHMP2B", 1), 
#                                                           rep("sALS", 1), 
#                                                           rep("FTLD", 1),
#                                                           rep("VCP", 1),
#                                                           rep("FUS", 1),
#                                                           rep("SOD1", 1))),
# col = list(Condition = c("C9_Pat" =  "olivedrab1", 
#                          "FTLD_Pat" = "skyblue1", 
#                          "sALS_Pat" = "gold", 
#                          "VCP_Pat" = "orchid2", 
#                          "CHMP2B_Pat" = "rosybrown1", 
#                          "FUS_Pat" = "azure3", 
#                          "SOD1_Pat" = "firebrick1",
#                          show_legend = TRUE)))
# 
# col = c("darkgreen", "green", "white", "red", "darkred")
# 
# Heatmap(LFC,
#         name = "Pathprint Score",
#         col = col,
#         cluster_rows = TRUE,
#         cluster_columns = TRUE,
#         clustering_method_rows = "complete",
#         clustering_distance_rows = "euclidean",
#         show_row_names = FALSE,
#         show_row_dend = TRUE, 
#         # heatmap_legend_param = list(at = c(-4,-2,0,1,3)),
#         # labels = c("-4", "-2", "0", "1", "3"),
#         # top_annotation = Disease,
#         column_names_max_height = unit(50, "mm"))
