##Pearson Correlation Coefficient for Expression Matrix

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15")
C9.GeneExpression <- read.csv(file = "C9rankeduniqueresult.csv")

ENS <- C9.GeneExpression[,15]

C9.disease <- C9.GeneExpression[,52:59]
C9.disease <- cbind(ENS, C9.disease)

rownames(C9.disease) <- C9.disease[,1]
C9.disease[,1] <- NULL
tC9.disease <- t(C9.disease)

C9.cor <- corr.test(tC9.disease, y= NULL, use = "pairwise", method = "pearson", adjust = "fdr")

DEG <- C9.cor[c("ACTN1", "ANXA1", "BBIP1", "BGN", "BPTF", "CDH11", "CREG1", "CSRP1", "CST3", "DCN", "GBAS", "JAG1",
                    "KCTD12", "KPNA6", "MPHOSPH9", "MXI1", "NDUFS5 /// RPL10", "NKTR", "NUTF2 /// NUTF2P4", "LOC101927673 /// OTUB1", 
                    "PFDN1", "PLEKHB1", "PLOD2", "POGZ", "PPP1R7", "PPP2CA", "PPP2CB", "PRPF3", "PTPN13", "RAB40B", "RPL35A", "RPL37", 
                    "SCN1B", "SERBP1", "SF3B1", "SPARC", "STMN1", "TARDBP", "TCF4", "TUBB4B", "TUG1", "VPS13B", "ZFP36", "ZFYVE26", "ZNF518A"),]
DEG <- t(DEG)

DEG.2 <- DEG[c("ACTN1", "ANXA1", "BBIP1", "BGN", "BPTF", "CDH11", "CREG1", "CSRP1", "CST3", "DCN", "GBAS", "JAG1",
               "KCTD12", "KPNA6", "MPHOSPH9", "MXI1", "NDUFS5 /// RPL10", "NKTR", "NUTF2 /// NUTF2P4", "LOC101927673 /// OTUB1", 
               "PFDN1", "PLEKHB1", "PLOD2", "POGZ", "PPP1R7", "PPP2CA", "PPP2CB", "PRPF3", "PTPN13", "RAB40B", "RPL35A", "RPL37", 
               "SCN1B", "SERBP1", "SF3B1", "SPARC", "STMN1", "TARDBP", "TCF4", "TUBB4B", "TUG1", "VPS13B", "ZFP36", "ZFYVE26", "ZNF518A"),]

DEG.2 <- as.data.frame(as.table(DEG.2))
DEG.2 <- subset(DEG.2, subset=(Freq !="1"))
DEG.2 <- DEG.2[!duplicated(DEG.2[,3]),] #check the remaining number is as expected

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/CorrelationNetwork/C9orf72")
write.csv(DEG.2, "SeedGenes.csv")
