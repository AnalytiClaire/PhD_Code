
setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/MedianGenes/")

## Load all datasets
C9 <- read.csv("C9rankeduniqueresult.csv")
C9 <- C9[order(C9$P.Value),]
CH <- read.csv("CHrankeduniqueresult.csv")
CH <- CH[order(CH$P.Value),]
sals <- read.csv("sALSrankeduniqueresult.csv")
sals <- sals[order(sals$P.Value),]
ftld <- read.csv("FTLDrankeduniqueresult.csv")
ftld <- ftld[order(ftld$P.Value),]
vcp <- read.csv("VCPrankeduniqueresult.csv")
vcp <- vcp[order(vcp$P.Value),]

pet <- read.csv("PET_HGNCENSrankeduniqueresult.csv")
rav <- read.csv("RAV_ensrankeduniqueresult.csv")


#FIND GENE NAMES FOR RAV DATA
ravgenes <- as.list(rav$Row.names) #wouldn't work as vector for some reason
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back_rav <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=ravgenes,  mart=mart)

#FIND GENE NAMES FOR PET DATA
genes <- as.vector(pet$Row.names)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back_pet <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=ravgenes,  mart=mart)

#Combine
pet <- merge(pet, mart_back_pet, by.x = "Row.names", by.y = "ensembl_gene_id")
pet <- pet[order(pet$P.Value),]
rav <- merge(rav, mart_back_rav, by.x = "Row.names", by.y = "ensembl_gene_id")
rav <- rav[order(rav$P.Value),]

## extract gene lists
c9_gene <- C9$Gene.Symbol
ch_gene <- CH$Gene.Symbol
sals_gene <- sals$Gene.Symbol
ftld_gene <- ftld$Gene.Symbol
vcp_gene <- vcp$Gene.Symbol
pet_gene <- pet$hgnc_symbol
rav_gene <- rav$hgnc_symbol



# num_overlap <- matrix(data=NA)
List <- list()

for (i in 1:6500){
  C9_int <- c9_gene[1:i]
  CH_int <- ch_gene[1:i]
  sals_int <- sals_gene[1:i]
  ftld_int <- ftld_gene[1:i]
  vcp_int <- vcp_gene[1:i]
  pet_int <- pet_gene[1:i]
  rav_int <- rav_gene[1:i]
  List[[i]] <- Reduce(intersect, list(C9_int, CH_int, sals_int, ftld_int, vcp_int, pet_int, rav_int))
}

# intersect = do.call(rbind, List)
# List <- as.matrix(List)
# List <- t(List)
# output <- matrix(unlist(List), nrow = 6000, byrow = TRUE)
output_6500 <- plyr::ldply(List, rbind)

write.csv(output_6500, "intersectnomedian_6000.csv")
write.csv(List, "list.csv",quote = FALSE, row.names = FALSE)


List[6500]

turnwhich(List == "165")
List[5877]
  


##### After changing the microarray gene names to only have one gene name (the first) ######
setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/MedianGenes/")

## Load all datasets
C9 <- read.csv("C9rankeduniqueresult.csv")
C9 <- C9[order(C9$P.Value),]
CH <- read.csv("CHrankeduniqueresult.csv")
CH <- CH[order(CH$P.Value),]
sals <- read.csv("sALSrankeduniqueresult.csv")
sals <- sals[order(sals$P.Value),]
ftld <- read.csv("FTLDrankeduniqueresult.csv")
ftld <- ftld[order(ftld),]
vcp <- read.csv("VCPrankeduniqueresult.csv")
vcp <- vcp[order(vcp$P.Value),]

pet <- read.csv("PET_HGNCENSrankeduniqueresult.csv")
rav <- read.csv("RAV_ensrankeduniqueresult.csv")


#FIND GENE NAMES FOR RAV DATA
ravgenes <- as.list(rav$Row.names) #wouldn't work as vector for some reason
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back_rav <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=ravgenes,  mart=mart)

#FIND GENE NAMES FOR PET DATA
genes <- as.vector(pet$Row.names)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back_pet <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=ravgenes,  mart=mart)

#Combine
pet <- merge(pet, mart_back_pet, by.x = "Row.names", by.y = "ensembl_gene_id")
pet <- pet[order(pet$P.Value),]
rav <- merge(rav, mart_back_rav, by.x = "Row.names", by.y = "ensembl_gene_id")
rav <- rav[order(rav$P.Value),]

## extract gene lists
c9_gene <- C9$Gene.Symbol
ch_gene <- CH$Gene.Symbol
sals_gene <- sals$Gene.Symbol
ftld_gene <- ftld$Gene.Symbol
vcp_gene <- vcp$Gene.Symbol
pet_gene <- pet$hgnc_symbol
rav_gene <- rav$hgnc_symbol



# num_overlap <- matrix(data=NA)
List <- list()

for (i in 1:6500){
  C9_int <- c9_gene[1:i]
  CH_int <- ch_gene[1:i]
  sals_int <- sals_gene[1:i]
  ftld_int <- ftld_gene[1:i]
  vcp_int <- vcp_gene[1:i]
  pet_int <- pet_gene[1:i]
  rav_int <- rav_gene[1:i]
  List[[i]] <- Reduce(intersect, list(C9_int, CH_int, sals_int, ftld_int, vcp_int, pet_int, rav_int))
}

# intersect = do.call(rbind, List)
# List <- as.matrix(List)
# List <- t(List)
# output <- matrix(unlist(List), nrow = 6000, byrow = TRUE)
output_6500 <- plyr::ldply(List, rbind)

write.csv(output_6500, "intersectnomedian_6000.csv")
write.csv(List, "list.csv",quote = FALSE, row.names = FALSE)


List[6500]

turnwhich(List == "165")
List[5877]
