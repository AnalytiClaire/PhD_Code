setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")

PPI <- read.table("iref14_Human_UP_noDup_table_nodash.txt", header = T)
braingenes <- read.csv("Zhang_BrainCelltype_Markers_braingenes.csv", header = T)
genelist <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/C9only_genenames.txt")
genelist <- genelist$V1

# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
# attributes <- listAttributes(mart)
# mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=genelist,  mart=mart)
# 
# genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
# setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/GRN")
# write.csv(genelist_Uniprot, "martback.csv", row.names = F)
# ###### IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING#####
# 
# 
# mart_table <- read.csv("martback.csv", header = T) #A table with the uniprot codes for the 285 DEGs
# uniprot_gene <- mart_table$uniprotswissprot
# 
# DEG_PPI <- subset(PPI, PPI$V1 %in% uniprot_gene | PPI$V2 %in% uniprot_gene)
# rownames(DEG_PPI) <- 1:nrow(DEG_PPI)
# 
# write.csv(DEG_PPI, "DEG_PPI.csv", row.names = F)

C9_uniprot <- "Q96LT7" #C9orf72
VCP_uniprot <- "P55072" #VCP
GRN_uniprot <- "P28799" #GRN

#First Neighbour
C9_PPI <- subset(PPI, PPI$V1 %in% C9_uniprot | PPI$V2 %in% C9_uniprot)
rownames(C9_PPI) <- 1:nrow(C9_PPI)

VCP_PPI <- subset(PPI, PPI$V1 %in% VCP_uniprot | PPI$V2 %in% VCP_uniprot)
rownames(VCP_PPI) <- 1:nrow(VCP_PPI)

GRN_PPI <- subset(PPI, PPI$V1 %in% GRN_uniprot | PPI$V2 %in% GRN_uniprot)
rownames(GRN_PPI) <- 1:nrow(GRN_PPI)

C91 <- as.matrix(C9_PPI$V1)
C92 <- as.matrix(C9_PPI$V2)
C9genes <- rbind(C91, C92)
C9genes <- unique(C9genes)

VCP1 <- as.matrix(VCP_PPI$V1)
VCP2 <- as.matrix(VCP_PPI$V2)
VCPgenes <- rbind(VCP1, VCP2)
VCPgenes <- unique(VCPgenes)

GRN1 <- as.matrix(GRN_PPI$V1)
GRN2 <- as.matrix(GRN_PPI$V2)
GRNgenes <- rbind(GRN1, GRN2)
GRNgenes <- unique(GRNgenes)

#Second Neighbour
C9_PPI_2 <- subset(PPI, PPI$V1 %in% C9genes | PPI$V2 %in% C9genes)
rownames(C9_PPI_2) <- 1:nrow(C9_PPI_2)

VCP_PPI_2 <- subset(PPI, PPI$V1 %in% VCPgenes | PPI$V2 %in% VCPgenes)
rownames(VCP_PPI_2) <- 1:nrow(VCP_PPI_2)

GRN_PPI_2 <- subset(PPI, PPI$V1 %in% GRNgenes | PPI$V2 %in% GRNgenes)
rownames(GRN_PPI_2) <- 1:nrow(GRN_PPI_2)
