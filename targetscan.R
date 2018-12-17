
PPIgene <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes.txt")
TS <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Predicted_Targets_Info.default_predictions_1miRNA.csv")


TSPPI <- subset(TS, TS$Gene.Symbol %in% PPIgene)

PPImiRNA <- unique(TS$miR.Family)
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/miRNA")
write.csv(TSPPI, "TargetScan_PPIhits.csv")
write.table(PPImiRNA, "uniquehits.txt", col.names = F, row.names = F, quote = F)
