# library(biomaRt)
# 
# setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/")
# 
# iref14 <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/iref14_Human_UP_noDup_table_nodash.txt", header = T)
# braingenes <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv", header = T)
# 
# genelist <- readLines("~/Documents/PhD/Parkinsons/ParkinsonsDiseaseMalacards.txt")
# 
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
# attributes <- listAttributes(mart)
# mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=genelist,  mart=mart)
# 
# setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/RPTS/")
# genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
# swiss <- subset(genelist_Uniprot, genelist_Uniprot$hgnc_symbol %in% genelist)
# write.csv(swiss, "martback.csv", row.names = F)
# 
# ###### IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING#####
# 
# mart_table <- read.csv("martback.csv", header = F)
# uniprot_gene <- mart_table$V2
# 
# PPI_All <- subset(iref14, iref14$V1 %in% uniprot_gene | iref14$V2 %in% uniprot_gene)
# write.csv(PPI_All, "PDmalacardsgene_PPI.csv", row.names = F, quote = F)
# 
# ## Convert Uniprot ID to HGNC symbol. Biomart jumbles output so
# ## Go To https://biodbnet-abcc.ncifcrf.gov/db/db2db.php submit each column of names. Select Uniprot Accession for input and Gene Symbol for output.
# # Select "NO" for remove duplicate input valies

setwd("~/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/RPTS/")
disgenePPI <- read.csv("PDmalacardsgene_PPI.csv")
genelist <- readLines("~/Documents/PhD/Parkinsons/ParkinsonsDiseaseMalacards.txt")

#make list of all genes in malacard network
PPI1 <- disgenePPI$Gene1
PPI1 <- as.character(unique(PPI1))
PPI2 <- disgenePPI$Gene2
PPI2 <- as.character(unique(PPI2))

PPI <- c(PPI1, PPI2)
PPI <- unique(PPI)

#Create list of genes that interact WITH disease genes, not including disease genes themselves
PPI <- subset(PPI, !(PPI %in% genelist))

### Load file with all genes ####
library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

sym.genes <- t(sym.genes)

allgenes <- sym.genes[!duplicated(sym.genes),]



### RUN RPT ####

test <- numeric()
#number of repetitions
m = 1000000
#Number of nug-malacard edges in malacards interaction network
nug <- 117
#Number of total nug genes
size <- 212

for (i in 1:m){
  sample <- sample(allgenes, size = size)
  test1 <- subset(disgenePPI, disgenePPI$Gene1 %in% sample)
  test2 <- subset(disgenePPI, disgenePPI$Gene2 %in% sample)
  test1gene <- as.character(test1$Gene1)
  test2gene <- as.character(test2$Gene2)
  testgene <- unique(c(test1gene, test2gene))
  test[[i]] <- length(testgene)
}

morethan <- which(test > nug)  # count number of times r is larger than test value
result <- sum((length(morethan)+1))/(m+1) # calculate P value
result
mean(test)
range(test)

par(oma=c(0,0,2,0))
hist(test,
     xlim = range(0:140),
     main = NULL, 
     xlab = "Number of Edges with Disease Genes")
abline(v = 117, col = "red", lwd = 2)
