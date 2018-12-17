library(biomaRt)

IPA <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/IPA/UpstreamAnalysis_nugget.csv")
GE <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/GeneXplain/NuggetUpstreamRegulators_genes.txt")
nug <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NuggetGenes.txt")
Disgenes <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")

library(Hmisc)
GE <- sapply(GE, toupper)
IPA <- IPA$Master.Regulator

intersect(GE, IPA)

genelist  <- as.character(c(GE, IPA))

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/")
iref14 <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/iref14_Human_UP_noDup_table_nodash.txt", header = T)
braingenes <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv", header = T)

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=genelist,  mart=mart)

genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
swiss <- subset(genelist_Uniprot, genelist_Uniprot$hgnc_symbol %in% genelist)
write.csv(swiss, "martback.csv", row.names = F)

###### IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING#####

disgenes <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/GeneXplain/")
mart_table <- read.csv("martback.csv", header = T)

uniprot_gene <- mart_table$uniprotswissprot
PPI_All <- subset(iref14, iref14$V1 %in% uniprot_gene & iref14$V2 %in% uniprot_gene)

write.csv(PPI_All, "ALL_PPI.csv", row.names = F, quote = F)
#symbol conversion https://biodbnet-abcc.ncifcrf.gov/db/db2db.php Uniprot Accession --> Gene Symbol

write.table(GE, "NuggetUpstreamRegulators_genes.txt", quote = F, col.names = F, row.names = F)
write.table(genelist, "AllUpstreamRegulators.txt", quote = F, col.names = F, row.names = F)



