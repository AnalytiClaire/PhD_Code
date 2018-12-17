setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
GRN <- read.csv("GRN_FTLDrankeduniqueresult.csv")
VCP <- read.csv("vcp_unique.csv")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")

DEGPPI <- read.table("DEG_PPI_Genes.txt")
DEGPPI <- DEGPPI$V1

mut_in_DEG <- subset(GRN, GRN$Gene.Symbol %in% DEGPPI)
mut_w_var <- data.frame(mut_in_DEG$Gene.Symbol, mut_in_DEG$adj.P.Val)
colnames(mut_w_var) <- c("GeneSymbol", "GRN_adjpval")
mut_w_var$GRNrank <- 1:length(mut_w_var$GeneSymbol)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/heatnode/")
write.csv(mut_w_var, "GRN_overlap_adjpval.csv", row.names = F)
