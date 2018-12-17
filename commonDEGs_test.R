setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
CH <- read.csv("CH_unique.csv")
CH <- CH[order(CH$P.Value),]
sals <- read.csv("sals_unique.csv")
sals <- sals[order(sals$P.Value),]
ftld <- read.csv("ftld_unique.csv")
ftld <- ftld[order(ftld$P.Value),]
vcp <- read.csv("vcp_unique.csv")
vcp <- vcp[order(vcp$P.Value),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
pet <- pet[!duplicated(pet$hgnc_symbol),]
rav <- read.csv("RAV_results_keepfiltering.csv")
rav <- rav[!duplicated(rav$hgnc_symbol),]

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

output_6500 <- plyr::ldply(List, rbind)

write.csv(output_6500, "intersectnomedian_6000.csv")
write.csv(List, "list.csv",quote = FALSE, row.names = FALSE)


List[6500]

turnwhich(List == "165")
List[5877]


x <- as.vector(List[6500])
cat(x, sep = "\n")

write.csv(List, "List.csv")

write.table(x, "6500_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
