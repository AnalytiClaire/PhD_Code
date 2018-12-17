setwd("/Users/clairegreen/Downloads")
load("DPD.Union.2017.Hs.gs")

setwd("/Users/clairegreen/Desktop/")
module <- read.table("NuggetGenes.txt")

#Load file with all genes
library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

sym.genes <- t(sym.genes)

allgenes <- sym.genes[!duplicated(sym.genes),]

library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)

martback<- getBM(attributes =c("hgnc_symbol", "entrezgene"), filters="hgnc_symbol", values=module,  mart=mart)

## Run Hyperpathway
# run script
pathwayEnrichment <- hyperPathway(
  genelist = martback$entrezgene,
  geneset = DPD.Hs.gs,
  Nchip = length(allgenes)
)
