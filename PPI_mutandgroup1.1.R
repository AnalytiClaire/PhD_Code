### PD Malacards ####
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")

PPI <- read.table("iref14_Human_UP_noDup_table_nodash.txt", header = T)

DEG_list <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/NUG_TDPgene_PPI.txt")

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=DEG_list,  mart=mart)
genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/")
write.csv(genelist_Uniprot, "martbacknug.csv", row.names = F)

# IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING #

mart_table <- read.csv("martback.csv", header = T) #A table with the uniprot codes for the DEGs
uniprot_gene <- mart_table$uniprotswissprot

DEG_PPI <- subset(PPI, PPI$V1 %in% uniprot_gene & PPI$V2 %in% uniprot_gene)
rownames(DEG_PPI) <- 1:nrow(DEG_PPI)

write.csv(DEG_PPI, "TDPNug_Disgenes_PPI.csv", row.names = F)

## Convert Uniprot ID to HGNC symbol. Biomart jumbles output so
## Go To https://biodbnet-abcc.ncifcrf.gov/db/db2db.php submit each column of names. Select Uniprot Accession for input and Gene Symbol for output.
# Select "NO" for remove duplicate input valies

DEG_PPI <- read.csv("TDPNug_Disgenes_PPI.csv")

DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene1 !="-")
DEG_PPI <- subset(DEG_PPI, DEG_PPI$Gene2 !="-")

write.csv(DEG_PPI, "FinalPDPPI.csv", row.names = F, quote = F)






















# for (i in 1:length(disgenes)){
#   mutgene <- disgenes[i]
#   mutgene_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=mutgene,  mart=mart)
#   martplusgene <- rbind(mart_table, mutgene_back[1,])
#   uniprot_gene <- mart_table$uniprotswissprot
#   PPI2 <- subset(iref14, iref14$V1 %in% uniprot_gene & iref14$V2 %in% uniprot_gene)
#   mutgeneprot <- mutgene_back$uniprotswissprot
#   interactions <- subset(PPI, PPI$V1 %in% mutgeneprot | PPI$V2 %in% mutgeneprot)
#   write.csv(interactions, file = paste(mutgene, "_PPI.csv", sep = ""), row.names = F, quote = F)
# }

mutgene <- disgenes
mutgene_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=mutgene,  mart=mart)
martplusgene <- rbind(mart_table, mutgene_back)
genelist_Uniprot <- subset(martplusgene, !(martplusgene$uniprotswissprot == ""))
genelist_Uniprot <- subset(genelist_Uniprot,!(duplicated(genelist_Uniprot$hgnc_symbol)))
uniprot_gene <- genelist_Uniprot$uniprotswissprot
PPI_All <- subset(iref14, iref14$V1 %in% uniprot_gene & iref14$V2 %in% uniprot_gene)
write.csv(PPI_All, "ALL_PPI.csv", row.names = F, quote = F)

#symbol conversion https://biodbnet-abcc.ncifcrf.gov/db/db2db.php

### Group1.1 only
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/PPI_mutgenes")
mart_table <- read.csv("martback.csv", header = T)
group1.1 <- mart_table$uniprotswissprot
PPI <- subset(iref14, iref14$V1 %in% group1.1 & iref14$V2 %in% group1.1)
write.csv(PPI, "Group1.1_PPI.csv", row.names = F, quote = F)



