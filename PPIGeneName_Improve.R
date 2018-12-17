#Trying to improve gene name coverage
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)

C9 <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/C9uniquegene_samples.csv")
names_C9 <- C9$Gene.Symbol
mart_back_C9 <- getBM(attributes =c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol", values=names_C9,  mart=mart)
mergeC9 <- merge(C9, mart_back_C9, by.x = "Gene.Symbol", by.y = "hgnc_symbol")


sALS <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/sALSuniquegene_samples.csv")
names_sALS <- sALS$Gene.Symbol
mart_back_sALS <- getBM(attributes =c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol", values=names_sALS,  mart=mart)
mergesALS <- merge(sALS, mart_back_sALS, by.x = "Gene.Symbol", by.y = "hgnc_symbol")


FTLD <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/FTLDuniquegene_samples.csv")
names_FTLD <- FTLD$Gene.Symbol
mart_back_FTLD <- getBM(attributes =c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol", values=names_FTLD,  mart=mart)
mergeFTLD <- merge(FTLD, mart_back_FTLD, by.x = "Gene.Symbol", by.y = "hgnc_symbol")


VCP <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/VCPuniquegene_samples.csv")
names_VCP <- VCP$Gene.Symbol
mart_back_VCP <- getBM(attributes =c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol", values=names_VCP,  mart=mart)
mergeVCP <- merge(VCP, mart_back_VCP, by.x = "Gene.Symbol", by.y = "hgnc_symbol")


PET <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/PET_results_keepfiltering.csv")
RAV <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/RAV_results_keepfiltering.csv")


#Convert DEG PPI nodes to ensembl IDs
DEG_PPI <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/OLDPPI/DEG_PPI_Genes_nofib.txt")
mart_DEGPPI <- getBM(attributes =c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol", values=DEG_PPI,  mart=mart)

#Subset each dataset with these common names so they are all the same size
C9_PPI   <- subset(mergeC9,   mergeC9$ensembl_gene_id %in% mart_DEGPPI$ensembl_gene_id)
VCP_PPI  <- subset(mergeVCP,  mergeVCP$ensembl_gene_id %in% mart_DEGPPI$ensembl_gene_id)
FTLD_PPI <- subset(mergeFTLD, mergeFTLD$ensembl_gene_id %in% mart_DEGPPI$ensembl_gene_id)
sALS_PPI <- subset(mergesALS, mergesALS$ensembl_gene_id %in% mart_DEGPPI$ensembl_gene_id)
PET_PPI  <- subset(PET,       PET$Gene %in% mart_DEGPPI$ensembl_gene_id)
RAV_PPI  <- subset(RAV,       RAV$Gene %in% mart_DEGPPI$ensembl_gene_id)


#Find the gene names that all datasets have in common
DEG_com <- Reduce(intersect, list(C9_PPI$ensembl_gene_id, sALS_PPI$ensembl_gene_id, FTLD_PPI$ensembl_gene_id,
                                  VCP_PPI$ensembl_gene_id, PET_PPI$Gene, RAV_PPI$Gene))

mart_DEG_com <- getBM(attributes =c("hgnc_symbol", "ensembl_gene_id"), filters="ensembl_gene_id", values=DEG_com,  mart=mart)

