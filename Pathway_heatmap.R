#####
load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprinted.RData")
load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/non-TDP/consensus_pathway.RData")

#Select patient columns
pat_C9 <- pathprint_C9[,4:11]
pat_sals <- pathprint_sals[,4:10]
pat_ftld <- pathprint_ftld[,9:24]
pat_vcp <- pathprint_vcp[,4:10]
#Combine
pat_all_pp <- cbind(pat_C9, pat_sals,pat_ftld, pat_vcp)
#Run consensus
pat_consensus <- data.frame(consensusFingerprint(pat_all_pp, thresh))

#Select control columns
con_C9 <- pathprint_C9[,1:3]
con_sals <- pathprint_sals[,1:3]
con_ftld <- pathprint_ftld[,1:8]
con_vcp <- pathprint_vcp[,1:3]
#Combine
con_all_pp <- cbind(con_C9, con_sals,con_ftld, con_vcp)
#Run consensus
con_consensus <- data.frame(consensusFingerprint(con_all_pp, thresh))

all_pp <- cbind(pat_C9,pat_ftld,pat_sals,pat_vcp, con_C9,con_ftld,con_sals,con_vcp)
all_pp <- subset(all_pp, !(rowSums(all_pp)==0))

dys_pp <- all_pp[rownames(all_pp) %in% dyspathways_PP,]
dys_maj <- all_pp[rownames(all_pp) %in% dyspathways_MAJ,]

library(gplots)
library(heatmap3)
library(ggplot2)

hmcols<- colorRampPalette(c("green","white", "red"))(3)
col<-colByValue(all_pp,col = c("green", "white", "red"), breaks = c(-1,0,1))
col<-colByValue(dys_pp,col=colorRampPalette(c('green',
                                            'white','red'))(3))
col = c("green", "white", "red")

#####
# Condition=c("yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#             "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#             "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#             "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#             "yellow","yellow","yellow","yellow","yellow","yellow","yellow", 
#             "magenta","magenta","magenta","magenta","magenta","magenta","magenta",
#             "magenta","magenta","magenta","magenta","magenta","magenta","magenta",
#             "magenta","magenta","magenta")
# colsidecolors <- cbind(Condition=c("yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#                                         "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#                                         "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#                                         "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#                                         "yellow","yellow","yellow","yellow","yellow","yellow","yellow", 
#                                         "magenta","magenta","magenta","magenta","magenta","magenta","magenta",
#                                         "magenta","magenta","magenta","magenta","magenta","magenta","magenta",
#                                         "magenta","magenta","magenta"),
#                             Disease=c("dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4",
#                                       "dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4",
#                                       "cyan","cyan", "cyan","cyan","cyan","cyan","cyan","cyan",
#                                       "cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan",
#                                       "cornflowerblue","cornflowerblue","cornflowerblue",
#                                       "cornflowerblue","cornflowerblue","cornflowerblue",
#                                       "cornflowerblue","blue","blue","blue","blue","blue","blue","blue",
#                                       "dodgerblue4","dodgerblue4","dodgerblue4",
#                                       "cyan","cyan", "cyan","cyan","cyan","cyan","cyan","cyan",
#                                       "cornflowerblue","cornflowerblue","cornflowerblue",
#                                       "blue","blue","blue"))
#                             # Tissue=c("magenta", "magenta", "purple", "purple4", "purple", "magenta"))
# 
# considecolors <- c("dodgerblue4","dodgerblue4","dodgerblue4","blue","blue","blue",
#                            "cyan","cyan", "cyan","cyan","cyan","cyan","cyan","cyan",
#                            "cornflowerblue","cornflowerblue","cornflowerblue")
# 
# # Condition=c("yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
# #             "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
# #             "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
# #             "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
# #             "yellow","yellow","yellow","yellow","yellow","yellow","yellow", 
# #             "magenta","magenta","magenta","magenta","magenta","magenta","magenta",
# #             "magenta","magenta","magenta","magenta","magenta","magenta","magenta",
#             "magenta","magenta","magenta")

Condition = c(rep("yellow", 38), rep("magenta", 17))

heatmap3(all_pp, 
         col = col,
         Rowv = TRUE,
         Colv = TRUE,
         distfun = dist,
         hclustfun = hclust,
         scale = "row",
         labCol = colnames(dys_pp),
         ColSideColors = colsidecolors,
         ColSideWidth = 1,
         cexCol = 0.4,
         cexRow = 0.5,
         xlab = "Samples",
         ylab = "Pathways",
         legendfun=function()
           showLegend(legend=c("-1","0","1"),col=c("green","white","red"),cex=1.5))

heatmap.2(dys_pp,
          Rowv = TRUE,
          Colv = NA,
          scale = "none",
          distfun = dist,
          hclustfun = hclust,
          dendrogram = "row",
          col=col,
          breaks = c(-1.5,-0.5,0.5,1),
          ColSideColors = Condition,
          key = TRUE,
          cexCol = 1.2,
          cexRow = 1,
          srtCol = 45,
          trace = "none",
          margins = c(10,20))

# library(reshape2)
# samples <- colnames(dys_pp)
# pathways <- rownames(dys_pp)
# colors <- c("green", "white", "red")
# ggplot(melt(cbind(sample=rownames(dys_pp), dys_pp)), aes(x = pathways, y = sample, fill = factor(value)), + 
#          geom_tile() + 
#          scale_fill_manual(values=colors))





######### COMPLEX HEATMAP ##########
library(ComplexHeatmap)

load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/allpathprint.RData")

#Select patient columns
pat_C9 <- pathprint_C9[,4:11]
pat_sals <- pathprint_sals[,4:10]
pat_ftld <- pathprint_ftld[,9:24]
pat_vcp <- pathprint_vcp[,4:10]
# pat_ch <- pathprint_chmp2b[,7:9]
pat_sod1 <- pathprint_sod1[,8:10]
pat_fus <- pathprint_fus[,4:6]
pat_CBFTLD <- pathprint_cbftld[,8:16]
#Combine
pat_all_pp <- cbind(pat_C9, pat_sals,pat_ftld, pat_vcp, pat_fus, pat_sod1, pat_CBFTLD)
#Run consensus
pat_consensus <- data.frame(consensusFingerprint(pat_all_pp, thresh))

#Select control columns
con_C9 <- pathprint_C9[,1:3]
con_sals <- pathprint_sals[,1:3]
con_ftld <- pathprint_ftld[,1:8]
con_vcp <- pathprint_vcp[,1:3]
# con_ch <- pathprint_chmp2b[,1:6]
con_sod1 <- pathprint_sod1[,1:7]
con_fus <- pathprint_fus[,1:3]
con_CBFTLD <- pathprint_cbftld[,1:7]
#Combine
con_all_pp <- cbind(con_C9, con_sals,con_ftld, con_vcp, con_fus, con_sod1, con_CBFTLD)
#Run consensus
con_consensus <- data.frame(consensusFingerprint(con_all_pp, thresh))

all_pp <- cbind(pat_C9,pat_ftld,pat_sals,pat_vcp, pat_fus, pat_sod1,pat_CBFTLD,
                con_C9,con_ftld,con_sals,con_vcp, con_fus, con_sod1, con_CBFTLD)
all_pp <- subset(all_pp, !(rowSums(all_pp)==0))

dyspathways_PP <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Consensus_fingerprint/PP_Filtered.txt", header = F)
dyspathways_PP <- dyspathways_PP$V1
dys_pp <- all_pp[rownames(all_pp) %in% dyspathways_PP,]

dys_pat <- pat_all_pp[rownames(pat_all_pp) %in% dyspathways_PP,]


col = c("cyan", "white", "red")

Condition = HeatmapAnnotation(df = data.frame(type1 = c(rep("Patient", 53), rep("Control", 34))), 
                              col = list(type1 = c("Patient" =  "red", "Control" = "blue")), show_legend = TRUE)
Disease = HeatmapAnnotation(df = data.frame(Condition = c(rep("C9_Pat", 8), 
                                                          rep("FTLD_Pat", 16), 
                                                          rep("sALS_Pat", 7), 
                                                          rep("VCP_Pat", 7),
                                                          rep("FUS_Pat", 3),
                                                          rep("SOD1_Pat", 3),
                                                          rep("CBFTLD_Pat",9),
                                                          rep("C9_Con", 3), 
                                                          rep("FTLD_Con", 8), 
                                                          rep("sALS_Con", 3), 
                                                          rep("VCP_Con", 3),
                                                          rep("FUS_Con", 3),
                                                          rep("SOD1_Con", 7),
                                                          rep("CBFTLD_Con",7))), 
                              col = list(Condition = c("C9_Pat" =  "olivedrab1", 
                                                       "C9_Con" = "olivedrab3", 
                                                       "FTLD_Pat" = "skyblue1", 
                                                       "FTLD_Con" = "skyblue3",
                                                       "sALS_Pat" = "gold", 
                                                       "sALS_Con" = "gold3",
                                                       "VCP_Pat" = "orchid2", 
                                                       "VCP_Con" = "orchid3",
                                                       "CHMP2B_Pat" = "rosybrown1", 
                                                       "CHMP2B_Con" = "rosybrown3",
                                                       "FUS_Pat" = "azure3", 
                                                       "FUS_Con" = "azure4",
                                                       "SOD1_Pat" = "firebrick1", 
                                                       "SOD1_Con" = "firebrick3",
                                                       "CBFTLD_Pat" = "royalblue1",
                                                       "CBFTLD_Con"= "navy")), show_legend = TRUE)

Heatmap(dys_pp,
        name = "Pathprint Score",
        col = col,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        clustering_method_rows = "complete",
        clustering_distance_rows = "euclidean",
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = TRUE,
        top_annotation = Disease, 
        column_names_max_height = unit(60, "mm"))


###########################################################################################################################
library(ComplexHeatmap)

load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/allpathprint.RData")

dyspathways_PP <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Consensus_fingerprint/PP_Filtered.txt", header = F)
dyspathways_PP <- dyspathways_PP$V1
# dyspathways_PP <- overlap

#Select patient columns
pat_C9 <- pathprint_C9[,4:11]
pat_C9_dys <- pat_C9[rownames(pat_C9) %in% dyspathways_PP,]
pat_sals <- pathprint_sals[,4:10]
pat_sals_dys <- pat_sals[rownames(pat_sals) %in% dyspathways_PP,]
pat_ftld <- pathprint_ftld[,9:24]
pat_ftld_dys <- pat_ftld[rownames(pat_ftld) %in% dyspathways_PP,]
pat_vcp <- pathprint_vcp[,4:10]
pat_vcp_dys <- pat_vcp[rownames(pat_vcp) %in% dyspathways_PP,]
# pat_ch <- pathprint_chmp2b[,7:9]
pat_sod1 <- pathprint_sod1[,8:10]
pat_sod1_dys <- pat_sod1[rownames(pat_sod1) %in% dyspathways_PP,]
pat_fus <- pathprint_fus[,4:6]
pat_fus_dys <- pat_fus[rownames(pat_fus) %in% dyspathways_PP,]
pat_CBFTLD <- pathprint_cbftld[,8:16]
pat_CBFTLD_dys <- pat_CBFTLD[rownames(pat_CBFTLD) %in% dyspathways_PP,]


#Select control columns
con_C9 <- pathprint_C9[,1:3]
con_C9_dys<- con_C9[rownames(con_C9) %in% dyspathways_PP,]
con_sals <- pathprint_sals[,1:3]
con_sals_dys<- con_sals[rownames(con_sals) %in% dyspathways_PP,]
con_ftld <- pathprint_ftld[,1:8]
con_ftld_dys<- con_ftld[rownames(con_ftld) %in% dyspathways_PP,]
con_vcp <- pathprint_vcp[,1:3]
con_vcp_dys<- con_vcp[rownames(con_vcp) %in% dyspathways_PP,]
# con_ch <- pathprint_chmp2b[,1:6]
con_sod1 <- pathprint_sod1[,1:7]
con_sod1_dys<- con_sod1[rownames(con_sod1) %in% dyspathways_PP,]
con_fus <- pathprint_fus[,1:3]
con_fus_dys<- con_fus[rownames(con_fus) %in% dyspathways_PP,]
con_CBFTLD <- pathprint_cbftld[,1:7]
con_CBFTLD_dys<- con_CBFTLD[rownames(con_CBFTLD) %in% dyspathways_PP,]




C9_pat_mean <- apply(pat_C9_dys, 1, mean)
C9_con_mean <- apply(con_C9_dys, 1, mean)
C9diff <- C9_pat_mean - C9_con_mean

sals_pat_mean <- apply(pat_sals_dys, 1, mean)
sals_con_mean <- apply(con_sals_dys, 1, mean)
salsdiff <- sals_pat_mean - sals_con_mean

ftld_pat_mean <- apply(pat_ftld_dys, 1, mean)
ftld_con_mean <- apply(con_ftld_dys, 1, mean)
ftlddiff <- ftld_pat_mean - ftld_con_mean

vcp_pat_mean <- apply(pat_vcp_dys, 1, mean)
vcp_con_mean <- apply(con_vcp_dys, 1, mean)
vcpdiff <- vcp_pat_mean - vcp_con_mean

sod1_pat_mean <- apply(pat_sod1_dys, 1, mean)
sod1_con_mean <- apply(con_sod1_dys, 1, mean)
sod1diff <- sod1_pat_mean - sod1_con_mean

fus_pat_mean <- apply(pat_fus_dys, 1, mean)
fus_con_mean <- apply(con_fus_dys, 1, mean)
fusdiff <- fus_pat_mean - fus_con_mean

cbftld_pat_mean <- apply(pat_CBFTLD_dys, 1, mean)
cbftld_con_mean <- apply(con_CBFTLD_dys, 1, mean)
cbftlddiff <- cbftld_pat_mean - cbftld_con_mean


PP_diffs <- data.frame(row.names = rownames(pat_C9_dys),
                       C9 = C9diff,
                       sals=salsdiff,
                       FTLD=ftlddiff,
                       VCP=vcpdiff,
                       SOD1=sod1diff,
                       FUS=fusdiff,
                       CBFTLD=cbftlddiff)



colsidecolorsLFC <- cbind(Disease=c(rep("skyblue4",2),
                                    "skyblue3",
                                    "skyblue1",
                                    rep("skyblue4",2),
                                    "skyblue3"),
                          Tissue=c(rep("gray25",2),
                                   "gray60",
                                   "gray80",
                                   rep("gray25",2),
                                   "gray95"),
                          Pathology=c(rep("magenta4",4),
                                      rep("magenta2", 3)))

hmcols<- colorRampPalette(c("green4","green","white", "red","red4"))(256)
heatmap3(PP_diffs, 
         col = hmcols,
         Rowv = TRUE,
         Colv = TRUE,
         distfun = dist,
         # hclustfun = hclust,
         scale = "row",
         labCol = colnames(PP_diffs),
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




