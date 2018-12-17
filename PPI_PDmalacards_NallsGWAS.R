library(biomaRt)

### PD Malacards ####
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")

PPI <- read.table("iref14_Human_UP_noDup_table_nodash.txt", header = T)
braingenes <- read.csv("Zhang_BrainCelltype_Markers_braingenes.csv", header = T)

DEG_list <- readLines("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PDNug_PDmalacards_PPI.txt")

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=DEG_list,  mart=mart)
genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/")
write.csv(genelist_Uniprot, "martbacknug.csv", row.names = F)

# IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING #

mart_table <- read.csv("martbacknug.csv", header = T) #A table with the uniprot codes for the DEGs
uniprot_gene <- mart_table$uniprotswissprot

DEG_PPI <- subset(PPI, PPI$V1 %in% uniprot_gene & PPI$V2 %in% uniprot_gene)
rownames(DEG_PPI) <- 1:nrow(DEG_PPI)

write.csv(DEG_PPI, "PDNug_PDmalacards_PPI.csv", row.names = F)

## Convert Uniprot ID to HGNC symbol. Biomart jumbles output so
## Go To https://biodbnet-abcc.ncifcrf.gov/db/db2db.php submit each column of names. Select Uniprot Accession for input and Gene Symbol for output.
# Select "NO" for remove duplicate input valies

DEG_PPI <- read.csv("PDNug_PDmalacards_PPI.csv")

DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene1 !="-")
DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene2 !="-")

write.csv(DEG_PPI, "FinalPDPPI.csv", row.names = F, quote = F)


### Nalls GWAS ####
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")

PPI <- read.table("iref14_Human_UP_noDup_table_nodash.txt", header = T)
braingenes <- read.csv("Zhang_BrainCelltype_Markers_braingenes.csv", header = T)

DEG_list <- readLines("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/Nugget_Nalls.txt")

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=DEG_list,  mart=mart)
genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/")
write.csv(genelist_Uniprot, "martbacknug.csv", row.names = F)

# IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING #

mart_table <- read.csv("martbacknug.csv", header = T) #A table with the uniprot codes for the DEGs
uniprot_gene <- mart_table$uniprotswissprot

DEG_PPI <- subset(PPI, PPI$V1 %in% uniprot_gene & PPI$V2 %in% uniprot_gene)
rownames(DEG_PPI) <- 1:nrow(DEG_PPI)

write.csv(DEG_PPI, "0point3_nallsPPI.csv", row.names = F)

## Convert Uniprot ID to HGNC symbol. Biomart jumbles output so
## Go To https://biodbnet-abcc.ncifcrf.gov/db/db2db.php submit each column of names. Select Uniprot Accession for input and Gene Symbol for output.
# Select "NO" for remove duplicate input valies

DEG_PPI <- read.csv("DEG_PPI_famBlood.csv")

DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene1 !="-")
DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene2 !="-")

write.csv(DEG_PPI, "FinalPDPPI.csv", row.names = F, quote = F)