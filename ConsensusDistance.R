##ConsensusDistance##
library(pathprint)
data(list = c("chipframe", "genesets","pathprint.Hs.gs","platform.thresholds", "pluripotents.frame"))

load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/allpathprint2.RData")

thresh = 0.4
#Select patient columns
pat_C9 <- pathprint_C9[,4:11]
pat_sals <- pathprint_sals[,4:10]
pat_GRN <- pathprint_GRN
pat_sftld <- pathprint_sftld
pat_vcp <- pathprint_vcp[,4:10]







#Run consensus
C9_consensus <- data.frame(consensusFingerprint(pat_C9, thresh))
sALS_consensus <- data.frame(consensusFingerprint(pat_sals, thresh))
GRN_consensus <- data.frame(consensusFingerprint(pat_GRN, thresh))
sFTLD_consensus <- data.frame(consensusFingerprint(pat_sftld, thresh))
VCP_consensus <- data.frame(consensusFingerprint(pat_vcp, thresh))

allconsensus <- c(C9_consensus, sALS_consensus, GRN_consensus, sFTLD_consensus, VCP_consensus)

pathprints <- list(pat_C9, pat_sals, pat_GRN, pat_sftld, pat_vcp)

consensusdistance <- list()


consensusdistance  <- consensusDistance(sALS_consensus, pat_C9)
mean(consensusdistance$distance)


CD <- matrix(0,5,5)

for (i in 1:length(allconsensus)){
  for (j in 1:length(pathprints)){
    cd  <- consensusDistance(allconsensus[[i]], pathprints[[j]])
    CD[i,j] <- mean(cd$distance) 
}}

rownames(CD) <- c("C9orf72_Con", "sALS_Con", "GRN_Con", "sFTLD_Con", "VCP_Con")
colnames(CD) <- c("C9orf72_AvgDis", "sALS_AvgDis", "GRN_AvgDis", "sFTLD_AvgDis", "VCP_AvgDis")
CD <- as.data.frame(CD)

# C9 vs Sals
AB <- mean(c(CD[1,2], CD[2,1]))
AC <- mean(c(CD[1,3], CD[3,1]))
AD <- mean(c(CD[1,4], CD[4,1]))
AE <- mean(c(CD[1,5], CD[5,1]))
BC <- mean(c(CD[2,3], CD[3,2]))
BD <- mean(c(CD[2,4], CD[4,2]))
BE <- mean(c(CD[2,5], CD[5,2]))
C.D <- mean(c(CD[3,4], CD[4,3]))
CE <- mean(c(CD[3,5], CD[5,3]))
DE <- mean(c(CD[4,5], CD[5,4]))

distances <- c(AB, AC, AD, AE, BC, BD, BE, C.D, CE, DE)
C1 <- c(rep("C9orf72", 4), rep("sALS", 3), rep("GRN", 2), "sFTLD")
C2 <- c("sALS", "GRN", "sFTLD", "VCP", "GRN", "sFTLD", "VCP","sFTLD", "VCP", "VCP")

Dist4Cyt <- data.frame(Dataset1 = C1,
                       Dataset2 = C2, 
                       Distance = distances)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Distance")
write.csv(Dist4Cyt, "consensusDistance_datasets_patonly.csv", row.names = F, quote = F)


#Calculate distance from the consensus
Distance<-consensusDistance(pathprint_C9, pathprint_sals)

#Plot histogram
# plot histograms
# par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
# C9Distance.hist<-hist(C9orf72Distance[,"distance"],
#                                    nclass = 50, xlim = c(0,1), main = "Distance from C9orf72 consensus")
# par(mar = c(7, 4, 4, 2))
# hist(C9orf72Distance[C9.LCM_pathprint, "distance"],
#      breaks = C9Distance.hist$breaks, xlim = c(0,1), 
#      main = "", xlab = "above: all GEO, below: curated C9orf72 samples")
