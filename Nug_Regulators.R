# Regulators

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/")

nug <- readLines("NuggetGenes.txt")
GE <- read.csv("GeneXplain/GEUpstreamRegulators.csv")
GE$Master.molecule.name <- toupper(GE$Master.molecule.name)
GE$Rank <- 1:nrow(GE)
IPA <- read.csv("IPA/UpstreamAnalysis_IPA.csv")
IPA$Rank <- 1:nrow(IPA)

GE_nug <- subset(GE, GE$Master.molecule.name %in% nug)
IPA_nug <- subset(IPA, IPA$Master.Regulator %in% nug)





#SNPNexus
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/ProjMine/")

nexus <- read.csv("NuggetSNPNexus1.0.csv")
info <- read.csv("SNP_info.csv")

info_nexus <- merge(nexus, info, by.x = "Start", by.y = "start")
meanpergene <- aggregate(info_nexus$p.val.K11., by=list(Category=info_nexus$gene_gemini), FUN=mean)

write.csv(info_nexus, "info_nexus.csv", row.names = F)
write.csv(meanpergene, "meanpergene_SNPnexus_K11Pval.csv", row.names = F)

