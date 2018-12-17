# ##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####
# 
library (pathprint)
data(list = c("chipframe", "genesets","pathprint.Hs.gs","platform.thresholds", "pluripotents.frame"))

options(stringsAsFactors = FALSE)
# 
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/NormalisedExpressionMatrices/") #set working directory to location of data
 
####C9_LCM ######
exp_C9.LCM <- read.csv ("C9eset.csv", header=TRUE, row.names = 1)
pathprint_C9 <- exprs2fingerprint(exp_C9.LCM, platform = "GPL570", species="human", progressBar=TRUE)

# ####sals_lcm###
exp_SALS.LCM <- read.csv ("sALSeset.csv", header=TRUE, row.names = 1)
pathprint_sals <- exprs2fingerprint (exp_SALS.LCM, platform = "GPL570", species="human", progressBar=TRUE)

# ####FTLD###
exp_FTLD <- read.csv ("FTLDeset.csv", header=TRUE, row.names = 1)
pathprint_ftld <- exprs2fingerprint(exp_FTLD, platform = "GPL571", species="human", progressBar=TRUE)

####VCP###
exp_VCP <- read.csv ("VCPeset.csv", header=TRUE, row.names = 1)
pathprint_vcp <- exprs2fingerprint (exp_VCP, platform = "GPL570", species="human", progressBar=T)
#

library(pathprint)

load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/allpathprint.RData")

load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint_biomart.RData")

thresh = 0.9
#Select patient columns
pat_C9 <- pathprint_C9[,4:11]
pat_sals <- pathprint_sals[,4:10]
pat_ftld <- pathprint_ftld[,9:24]
pat_vcp <- pathprint_vcp[,4:10]
#Combine
pat_all_pp <- cbind(pat_C9, pat_sals,pat_ftld, pat_vcp)

# 
# keep <- rowSums(pat_all_pp==1) >= ncol(pat_all_pp)/(100/80)
# pat_up<-pat_all_pp[keep,]
# keep <- rowSums(pat_all_pp==-1) >= ncol(pat_all_pp)/(100/80)
# pat_down<-pat_all_pp[keep,]


#Run consensus
pat_consensus <- data.frame(consensusFingerprint(pathprint_C9, thresh))

#Select control columns
con_C9 <- pathprint_C9[,1:3]
con_sals <- pathprint_sals[,1:3]
con_ftld <- pathprint_ftld[,1:8]
con_vcp <- pathprint_vcp[,1:3]
#Combine
con_all_pp <- cbind(con_C9, con_sals,con_ftld, con_vcp)

# keep <- rowSums(con_all_pp==1) >= ncol(con_all_pp)/(100/80)
# con_up<-con_all_pp[keep,]
# keep <- rowSums(con_all_pp==-1) >= ncol(con_all_pp)/(100/80)
# con_down<-con_all_pp[keep,]

#Run consensus
con_consensus <- data.frame(consensusFingerprint(con_all_pp, thresh))

# #Split up and downregulated pathways
# pat_up <- subset(pat_consensus, pat_consensus[,1] ==1)
# pat_down <- subset(pat_consensus, pat_consensus[,1]==-1)
# 
# con_up <- subset(con_consensus, con_consensus[,1]==1)
# con_down <- subset(con_consensus, con_consensus[,1]==-1)
# 
# #Extract pathway names
# up_pat_path <- rownames(pat_up)
# down_pat_path <- rownames(pat_down)
# up_con_path <- rownames(con_up)
# down_con_path <- rownames(con_down)
# 
# #Subset patient pathways that are not present in control
# x <- subset(pat_up, !(rownames(pat_up) %in% up_con_path))
# y <- subset(pat_down, !(rownames(pat_down) %in% down_con_path))
# 


##########
all_consensus <- merge(pat_consensus, con_consensus, by = 0)
diff_consensus_PP <- subset(all_consensus, !(all_consensus$consensusFingerprint.pathprint_C9..thresh. == all_consensus$consensusFingerprint.con_all_pp..thresh.))

down_path <- subset(diff_consensus_PP, diff_consensus_PP$consensusFingerprint.pathprint_C9..thresh. < diff_consensus_PP$consensusFingerprint.con_all_pp..thresh.)
down_pathways <- down_path$Row.names
up_path <- subset(diff_consensus_PP, diff_consensus_PP$consensusFingerprint.pathprint_C9..thresh. > diff_consensus_PP$consensusFingerprint.con_all_pp..thresh.)
up_pathways <- up_path$Row.names

dyspathways_PP <- c(down_pathways, up_pathways)




all_pp <- cbind(pat_C9,pat_ftld,pat_sals,pat_vcp, con_C9,con_ftld,con_sals,con_vcp)
# all_pp <- subset(all_pp, !(rowSums(all_pp)==0))

dys_pp <- all_pp[rownames(all_pp) %in% dyspathways_PP,]

patmean <- apply(pat_all_pp, 1, mean)
names(patmean) <- rownames(pat_all_pp)
conmean <- apply(con_all_pp, 1, mean)
names(conmean) <- rownames(con_all_pp)

mean.df <- data.frame(row.names = rownames(pat_all_pp),
                      Pat.Mean= patmean,
                      Con.Mean = conmean)
mean.df$difference <- (mean.df$Pat.Mean - mean.df$Con.Mean)

dys.mean<- subset(mean.df, rownames(mean.df) %in% dyspathways_PP)
dys.mean$difference <- (dys.mean$Pat.Mean - dys.mean$Con.Mean)
dys.mean <- dys.mean[order(dys.mean$difference, decreasing = TRUE),]
dys.mean <- dys.mean[ which( dys.mean$difference > 0.1 | dys.mean$difference < -0.1) , ]
list <- rownames(dys.mean)


cat(rownames(dys.mean), sep = '\n')
cat(dys.mean$Pat.Mean, sep = '\n')
cat(dys.mean$Con.Mean, sep = '\n')
cat(dys.mean$difference, sep = '\n')

dys.mean[14,] <- NULL

#####
# consensus_pat = subset(pat_all_pp, row.names(pat_all_pp) %in% dyspathways_PP)
# consensus_pat <- merge(consensus_pat, dys.mean, by =0)
# consensus_pat <- consensus_pat[order(consensus_pat$difference, decreasing = TRUE),]
# rownames(consensus_pat) <- consensus_pat$Row.names
# consensus_pat[,1] <- NULL
# consensus_pat <- t(consensus_pat[,1:38])
# 
# consensus_con = subset(con_all_pp, row.names(con_all_pp) %in% dyspathways_PP)
# consensus_con <- merge(consensus_con, dys.mean, by =0)
# consensus_con <- consensus_con[order(consensus_con$difference, decreasing = TRUE),]
# rownames(consensus_con) <- consensus_con$Row.names
# consensus_con[,1] <- NULL
# consensus_con <- t(consensus_con[,1:17])
# 
# boxplot(consensus_pat, las = 2,
#         col = c(rep("Red", 23), rep("Blue", 24)))
# boxplot(consensus_con, las = 2,
#         col = c(rep("Red", 23), rep("Blue", 24)))
# 
# setwd("~/Desktop/Desktop Folder/")
# write.csv(consensus_pat, "consensus_pat.csv")

#### USING CONSENSUS ####
# bar.dcp <- diff_consensus_PP
# bar.dcp <- merge(bar.dcp, dys.mean, by.x = "Row.names", by.y = 0)
# bar.dcp <- bar.dcp[order(bar.dcp$difference, decreasing = TRUE),]
# rownames(bar.dcp) <- bar.dcp$Row.names
# bar.dcp[,1] <- NULL
# bar.dcp <- t(bar.dcp[,1:2])
# 
# barpat <- t(as.data.frame(bar.dcp[1,]))
# barcon <- t(as.data.frame(bar.dcp[2,]))


#### USING MEAN ####
mean.df <- mean.df[order(mean.df$difference, decreasing = T),]
topmean <- mean.df[1:23,]
bottommean <- mean.df[610:633,]

pathways <- rbind(topmean,bottommean)
bar.mean <- merge(all_consensus, pathways, by.x = "Row.names", by.y = 0)



bar_PP <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Consensus_fingerprint/Barchart_pathprint.csv")
bar_PP_mean <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Consensus_fingerprint/Barchart_mean.csv")


eg_pat <- pat_all_pp["Complement and coagulation cascades (KEGG)",]
eg_con <- con_all_pp["Complement and coagulation cascades (KEGG)",]
eg_pat <- as.vector(eg_pat)
eg_con <- as.vector(eg_con)

eg <- as.vector(c(eg_pat,eg_con))

cols <- c("black", "red")

plot(eg, col = c(rep("black", 38), rep("red", 17)),
     xlab = "Sample", 
     ylab = "Pathprint Score",
     main = "Complement and coagulation cascades (KEGG)")
legend(50,1,legend = c("Patient","Control"),
       col = c("black", "red"),
       fill = NULL)

# plot(eg_pat, xlab = "Sample")
# points(eg_con, col = "red", pch = 8)





library(ggplot2)
ggplot(bar_PP_mean, aes(factor(Pathway, levels=unique(Pathway)), Value, fill = Condition)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# barplot(barpat, las = 2)


library(ComplexHeatmap)

col = c("green", "white", "red")

Condition = HeatmapAnnotation(df = data.frame(Condition = c(rep("Patient", 38), rep("Control", 17))), 
                              col = list(Condition = c("Patient" =  "red", "Control" = "blue")), show_legend = TRUE)
Disease = HeatmapAnnotation(df = data.frame(Condition = c(rep("C9_Pat", 8), rep("FTLD_Pat", 16), rep("sALS_Pat", 7), rep("VCP_Pat", 7),
                                                          rep("C9_Con", 3), rep("FTLD_Con", 8), rep("sALS_Con", 3), rep("VCP_Con", 3))), 
                            col = list(Condition = c("C9_Pat" =  "olivedrab1", "C9_Con" = "olivedrab4", 
                                                     "FTLD_Pat" = "skyblue1", "FTLD_Con" = "skyblue3",
                                                     "sALS_Pat" = "gold", "sALS_Con" = "goldenrod2",
                                                     "VCP_Pat" = "orchid2", "VCP_Con" = "orchid4")), show_legend = TRUE)

Heatmap(dys_pp,
        name = "Pathprint",
        col = col,
        cluster_rows = TRUE,
        cluster_columns = F,
        clustering_method_rows = "complete",
        clustering_distance_rows = "euclidean",
        show_row_names = FALSE,
        show_row_dend = TRUE,
        top_annotation = Condition, 
        column_names_max_height = unit(50, "mm"))


###################################################################################################################
################################## NON - TDP-43 ###################################################################
###################################################################################################################

library(pathprint)
library(affy)
library(tkWidgets)
# setwd("/Users/clairegreen/Documents/PhD/TDP-43/non-TDP-43 Data Sets/J-Y_FUS")
# #run program to choose .CEL files from directory
# celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
# #celfiles<-basename(celfiles)
# Data<-ReadAffy(filenames=celfiles) #read in files
# rmaEset<-rma(Data) #normalise using RMA
# FUS_norm<-exprs(rmaEset) #takes expression from normalised expression set
# 
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/FTD-U.brain/")
#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
#celfiles<-basename(celfiles)
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
CB_FTLD_norm<-exprs(rmaEset) #takes expression from normalised expression set
# 
# ####FUS###
# pathprint_fus <- exprs2fingerprint (FUS_norm, platform = "GPL570", species="human", progressBar=T)
# ####SOD1###
# pathprint_sod1 <- exprs2fingerprint (SOD1_norm, platform = "GPL570", species="human", progressBar=T)
####CB_FTLD###
pathprint_cbftld <- exprs2fingerprint (CB_FTLD_norm, platform = "GPL571", species="human", progressBar=T)
#####

load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/allpathprint.RData")
load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint_biomart.RData")

thresh = 0.4
#Select patient columns
pat_FUS <- pathprint_fus[,4:6]
pat_SOD1 <- pathprint_sod1[,8:10]
pat_CBFTLD <- pathprint_cbftld[,8:16] 

#Combine
pat_all_pp <- cbind(pat_FUS, pat_SOD1)

#Run consensus
pat_consensus <- data.frame(consensusFingerprint(pat_all_pp, thresh))

#Select control columns
con_FUS <- pathprint_fus[,1:3]
con_SOD1 <- pathprint_sod1[,1:7]
con_CBFTLD <- pathprint_cbftld[,1:7]
#Combine
con_all_pp <- cbind(con_FUS, con_SOD1)

#Run consensus
con_consensus <- data.frame(consensusFingerprint(con_all_pp, thresh))


##########
all_consensus <- merge(pat_consensus, con_consensus, by = 0)
diff_consensus_PP <- subset(all_consensus, !(all_consensus$consensusFingerprint.pat_all_pp..thresh. == all_consensus$consensusFingerprint.con_all_pp..thresh.))

down_path <- subset(diff_consensus_PP, diff_consensus_PP$consensusFingerprint.pat_all_pp..thresh. < diff_consensus_PP$consensusFingerprint.con_all_pp..thresh.)
down_pathways <- down_path$Row.names
up_path <- subset(diff_consensus_PP, diff_consensus_PP$consensusFingerprint.pat_all_pp..thresh. > diff_consensus_PP$consensusFingerprint.con_all_pp..thresh.)
up_pathways <- up_path$Row.names

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/non-TDP")
write.table(down_pathways, "down_pathways.txt", col.names = F, row.names = F, quote = F)
write.table(up_pathways, "up_pathways.txt", col.names = F, row.names = F, quote = F)
down_pathways
up_pathways

dyspathways_PP <- c(down_pathways, up_pathways)



patmean <- apply(pat_all_pp, 1, mean)
names(patmean) <- rownames(pat_all_pp)
conmean <- apply(con_all_pp, 1, mean)
names(conmean) <- rownames(con_all_pp)

mean.df <- data.frame(row.names = rownames(pat_all_pp),
                      Pat.Mean= patmean,
                      Con.Mean = conmean)
mean.df$difference <- (mean.df$Pat.Mean - mean.df$Con.Mean)

dys.mean.non<- subset(mean.df, rownames(mean.df) %in% dyspathways_PP)
dys.mean.non$difference <- (dys.mean.non$Pat.Mean - dys.mean.non$Con.Mean)
dys.mean.non <- dys.mean.non[order(dys.mean.non$difference, decreasing = TRUE),]
dys.mean.non <- dys.mean.non[ which( dys.mean.non$difference > 0.2 | dys.mean.non$difference < -0.2) , ]



cat(rownames(dys.mean), sep = '\n')
cat(dys.mean$Pat.Mean, sep = '\n')
cat(dys.mean$Con.Mean, sep = '\n')
cat(dys.mean$difference, sep = '\n')

###################################################################################################################
################################## FTLD Cerebellum ################################################################
###################################################################################################################

