### Arthritis PPI Network ####
library(biomaRt)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")

PPI <- read.table("iref14_Human_UP_noDup_table_nodash.txt", header = T)
braingenes <- read.csv("Zhang_BrainCelltype_Markers_braingenes.csv", header = T)

DEG_list <- readLines("~/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/allDEGs_SJO_SLE_ALO.txt")

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=DEG_list,  mart=mart)
genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PPI")
write.csv(genelist_Uniprot, "martback.csv", row.names = F)

# IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING #

mart_table <- read.csv("martback_Aug2018.csv", header = T) #A table with the uniprot codes for the DEGs
uniprot_gene <- mart_table$uniprotswissprot

DEG_PPI <- subset(PPI, PPI$V1 %in% uniprot_gene | PPI$V2 %in% uniprot_gene)
rownames(DEG_PPI) <- 1:nrow(DEG_PPI)

write.csv(DEG_PPI, "DEG_PPI_RAOA.csv", row.names = F)

## Convert Uniprot ID to HGNC symbol. Biomart jumbles output so
## Go To https://biodbnet-abcc.ncifcrf.gov/db/db2db.php submit each column of names. Select Uniprot Accession for input and Gene Symbol for output.
# Select "NO" for remove duplicate input valies

DEG_PPI <- read.csv("DEG_PPI_RAOA.csv")

DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene1 !="-")
DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene2 !="-")

write.csv(DEG_PPI, "FinalPDPPI.csv", row.names = F, quote = F)

## 6600 interactions ###

# AFTER DOWNLOADING NODE TABLE FROM CYTOSCAPE #

nodetable <- read.csv("FinalPDPPInode.csv")
PDgenes <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/ParkinsonsDiseaseMalacards.txt")
celltype <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv")
DEG <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/ALS_sfblood_ALLgenes.txt")

nodetable$PDMalacards <- nodetable$shared.name %in% PDgenes
nodetable$DEG <- nodetable$name %in% DEG
nodetable_celltype <- merge(celltype, nodetable, by.x = "Gene.symbol",  by.y = "shared.name", all = T)
nodetable_celltype <- subset(nodetable_celltype, !(nodetable_celltype$SUID == "NA"))

write.csv(nodetable_celltype, "ModifiedPDNodetable.csv", row.names = F)