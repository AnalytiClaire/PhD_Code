###Random permutations for Differential Expression##
options(scipen=999)

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")

#load the total number of unique genes for each data set
# Group1 <- read.csv("C9rankeduniqueresult.csv")
# Group2 <- read.csv("CHrankeduniqueresult.csv")
# Group3 <- read.csv("sALSrankeduniqueresult.csv")
# Group4 <- read.csv("FTLDrankeduniqueresult.csv")
# Group5 <- read.csv("VCPrankeduniqueresult.csv")
# Group6 <- read.csv("PETNOrankeduniqueresult.csv")
# Group7 <- read.csv("RAVrankeduniqueresult.csv")
# 

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
# CH <- read.csv("CH_unique.csv")
# CH <- CH[order(CH$P.Value),]
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
# ch_gene <- CH$Gene.Symbol
sals_gene <- sals$Gene.Symbol
ftld_gene <- ftld$Gene.Symbol
vcp_gene <- vcp$Gene.Symbol
pet_gene <- pet$hgnc_symbol
rav_gene <- rav$hgnc_symbol

#indicate the number of overlapping genes identified by DE analysis
test <- 328
# samplenum <- 6000

m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

## Gene numbers found by running "Foldchangeupdown.R"
for (j in 1:m){
  random1 <- sample (c9_gene, size=3788, replace=FALSE)
  # random2 <- sample (ch_gene, size=samplenum, replace=FALSE)
  random3 <- sample (sals_gene, size=5905, replace=FALSE)
  random4 <- sample (ftld_gene, size=4941, replace=FALSE)
  random5 <- sample (vcp_gene, size=8011, replace=FALSE)
  random6 <- sample (pet_gene, size=9259, replace=FALSE)
  random7 <- sample (rav_gene, size=8028, replace=FALSE)
  random <- Reduce(intersect, list(random1, random3, random4, random5, random6, random7))
  r[j] <- length(random)
}

test1 <- which(r > test)  # count number of times r is larger than test value
result <- (length(test1)/m) # calculate P value
result
mean(r)
range(r)
