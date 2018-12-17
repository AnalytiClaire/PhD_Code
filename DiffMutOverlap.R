##Consensus of genes generated from Coxpresdb

# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/CoxpresDB/CDB_Top500/")
# 
# GRN <- read.csv("GRN_500_coex_list.csv")
# OPTN <- read.csv("OPTN_500_coex_list.csv")
# SQSTM1 <- read.csv("SQSTM1_500_coex_list.csv")
# TARDBP <- read.csv("TARDBP_500_coex_list.csv")
# TBK1 <- read.csv("TBK1_500_coex_list.csv")
# UBQLN2 <- read.csv("UBQLN2_500_coex_list.csv")
# VCP <- read.csv("VCP_500_coex_list.csv")
# hnRNPA1 <- read.csv("hnRNPA1_500_coex_list.csv")
# 
# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/CoxpresDB/CDB_Top2000/")
# GRN <- read.csv("GRN_2000_coex_list.csv")
# OPTN <- read.csv("OPTN_2000_coex_list.csv")
# SQSTM1 <- read.csv("SQSTM1_2000_coex_list.csv")
# TARDBP <- read.csv("TARDBP_2000_coex_list.csv")
# TBK1 <- read.csv("TBK1_2000_coex_list.csv")
# UBQLN2 <- read.csv("UBQLN2_2000_coex_list.csv")
# VCP <- read.csv("VCP_2000_coex_list.csv")
# hnRNPA1 <- read.csv("hnRNPA1_2000_coex_list.csv")

# GRN_gene <- GRN$Gene
# OPTN_gene <- OPTN$Gene
# SQSTM1_gene <- SQSTM1$Gene
# TARDBP_gene <- TARDBP$Gene
# TBK1_gene <- TBK1$Gene
# UBQLN2_gene <- UBQLN2$Gene
# VCP_gene <- VCP$Gene
# hnRNPA1_gene <- hnRNPA1$Gene

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/CoxpresDB/CDB_Top500/")

GRN <- read.csv("analysis_GRN_500.csv")
OPTN <- read.csv("analysis_OPTN_500.csv")
SQSTM1 <- read.csv("analysis_SQSTM1_500.csv")
TARDBP <- read.csv("analysis_TARDBP_500.csv")
TBK1 <- read.csv("analysis_TBK1_500.csv")
UBQLN2 <- read.csv("analysis_UBQLN2_500.csv")
VCP <- read.csv("analysis_VCP_500.csv")
hnRNPA1 <- read.csv("analysis_hnRNPA1_500.csv")

GRN_gene <- GRN$GO.biological.process.complete
OPTN_gene <- OPTN$GO.biological.process.complete
SQSTM1_gene <- SQSTM1$GO.biological.process.complete
TARDBP_gene <- TARDBP$GO.biological.process.complete
TBK1_gene <- TBK1$GO.biological.process.complete
UBQLN2_gene <- UBQLN2$GO.biological.process.complete
VCP_gene <- VCP$GO.biological.process.complete
hnRNPA1_gene <- hnRNPA1$GO.biological.process.complete


overlap <- Reduce(intersect, list(GRN_gene, OPTN_gene, SQSTM1_gene, TARDBP_gene, TBK1_gene, UBQLN2_gene, VCP_gene, hnRNPA1_gene))
overlap
