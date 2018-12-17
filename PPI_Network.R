library(biomaRt)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/CurrentPPI")

PPI <- read.table("iref14_Human_UP_noDup_table_nodash.txt", header = T)
braingenes <- read.csv("Zhang_BrainCelltype_Markers_braingenes.csv", header = T)
DEG_list <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/CurrentStuff/Filtered_all_genes.txt")
DEG_list <- DEG_list$V1

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=DEG_list,  mart=mart)

genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/CurrentPPI/")
write.csv(genelist_Uniprot, "martback.csv", row.names = F)
###### IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING#####

mart_table <- read.csv("martback.csv", header = T) #A table with the uniprot codes for the DEGs
uniprot_gene <- mart_table$uniprotswissprot

DEG_PPI <- subset(PPI, PPI$V1 %in% uniprot_gene | PPI$V2 %in% uniprot_gene)
rownames(DEG_PPI) <- 1:nrow(DEG_PPI)

write.csv(DEG_PPI, "DEG_PPI.csv", row.names = F)


# V1 <- DEG_PPI$V1
# write.table(V1, "V1.txt", row.names = F, col.names = F, quote = F)
# V2 <- DEG_PPI$V2
# 
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
# mart_back_V1 <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="uniprotswissprot", values=V1,  mart=mart)
# 
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
# mart_back_V2 <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="uniprotswissprot", values=V2,  mart=mart)


### Biomart reorders the genes when converting so I used the online tool bioDBnet (https://biodbnet-abcc.ncifcrf.gov/db/db2db.php)
# using the UniProt accessions as input and selecting "Gene Symbol". MAKE SURE YOU DO NOT REMOVE DUPLICATE INPUT VALUES


DEG_PPI <- read.csv("DEG_PPI.csv")
DEG_PPI <- DEG_PPI[DEG_PPI$Gene1 !="-" & DEG_PPI$Gene2 !="-",]


# V1_table <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/GRN/V1_bioDBnet_db2db_170526051436_1018514005.txt", header = T)
# V2_table <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/GRN/V2_bioDBnet_db2db_170526051450_1047900660.txt", header = T)
# 
# DEG_Table <- cbind(V1_table, V2_table)
# DEG_Table <- DEG_Table[c(1,3,2,4)]
# 
# DEG_Table$OG_Gene <- DEG_Table$GeneSymbol %in% genelist 
# DEG_Table$OG_Gene <- DEG_Table$GeneSymbol.1 %in% genelist
# 
# CYT <- DEG_Table[,3:4]
# CYT$INT <- "PPI"
# write.csv(CYT, "CYT_genesymbol_GRN.csv", row.names = F, quote = F)



############# NOW GO LOAD INTO CYTOSCAPE ###############################




############# AFTER DOWNLOADING NODE TABLE FROM CYTOSCAPE ##############
nodetable <- read.csv("CYT_genesymbol_GRNdefaultnode.csv")
ALSOD <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/ALSOD.txt")
ALSOD <- ALSOD$V1

celltype <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv")
nodetable$OG_Gene <- nodetable$name %in% genelist
nodetable$ALSOD <- nodetable$shared.name %in% ALSOD
nodetable$DEG <- nodetable$name %in% DEG_list
nodetable_celltype <- merge(celltype, nodetable, by.x = "Gene.symbol",  by.y = "shared.name", all = T)
nodetable_celltype <- subset(nodetable_celltype, !(nodetable_celltype$SUID == "NA"))

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/GRN/")
write.csv(nodetable_celltype, "GRN_Node_Table.csv", row.names = F)




###### merging networks ######
C9node <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/C9/C9top100.txt", header = F)
C9node <- C9node$V1
VCPnode <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/VCP/VCPtop100.txt", header = F)
VCPnode <- VCPnode$V1
GRNnode <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/VCP/VCPtop100.txt", header = F)
GRNnode <- GRNnode$V1

mergenode <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/MergedNetwork_2defaultnode.csv")
mergenode$OG_Gene <- NULL
mergenode$C9 <- mergenode$name %in% C9node
mergenode$VCP <- mergenode$name %in% VCPnode

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/")
write.csv(mergenode, "mergenode.csv")
