##### Finding perfect cutoff with separate benchmarks

library(pathprint)
library(hgu133plus2.db)

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

for (i in 1:8000){
  C9_int <- c9_gene[1:i]
  CH_int <- ch_gene[1:i]
  sals_int <- sals_gene[1:i]
  ftld_int <- ftld_gene[1:i]
  vcp_int <- vcp_gene[1:i]
  pet_int <- pet_gene[1:i]
  rav_int <- rav_gene[1:i]
  List[[i]] <- Reduce(intersect, list(C9_int, CH_int, sals_int, ftld_int, vcp_int, pet_int, rav_int))
}

#Load file with all genes
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]
sym.genes <- t(sym.genes)
allgenes <- sym.genes[!duplicated(sym.genes),]


#Load Benchmark gene lists
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/")
W <- read.csv(file = "BenchmarkGenes.csv", na.strings = c("", "NA)"))
W <- as.list(W)
W<- lapply(W, function(x) x[!is.na(x)])


#Remove list elements with less than 5 genes (to aid calculations)
List_5 <- List[lengths(List) > 4]
#Leaves final 5087 elements (elements 1:2913 removed)

#Create new empty list
enrich_result <- list()

#Run for loop to calculate hyperPathway for different consensus thresholds
for (i in 1:length(List_5)){
  
  write.table(List_5[i], "benchmark_genelist.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  bgene <- read.table("benchmark_genelist.txt")
  bgene <- bgene$V1
  
  pathwayEnrichment <- hyperPathway(
    genelist = bgene,
    geneset = W,
    Nchip = length(allgenes))
  
  enrich_result[[i]] <- pathwayEnrichment
}

#Run for loop to extract each row into a separate list
Exac <- list()
Taylor <- list()
Pasterkamp <- list()
Subnetwork.28 <- list()
GeneCards.ALS <- list()
GeneCards.AD<- list()
GWASCentral.ALS <- list()
GWASCentral.AD <- list()
NeuroX.05 <- list()
NeuroX.GWS <- list()
ALSOD <- list()
Cirulli <- list()


#Collect enrichments for each list 
for (i in 1:length(enrich_result)){
  br <- data.frame(enrich_result[i],stringsAsFactors = F)
  br <- data.frame(t(br),stringsAsFactors = F)
  Exac[[i]] <- br$Exac
  Taylor[[i]] <- br$Taylor
  Pasterkamp[[i]] <- br$Pasterkamp
  Subnetwork.28[[i]] <- br$Subnetwork.28
  GeneCards.ALS[[i]] <- br$GeneCards.ALS
  GeneCards.AD[[i]]<- br$GeneCards.AD
  GWASCentral.ALS[[i]] <- br$GWASCentralALS
  GWASCentral.AD[[i]] <- br$GWASCentralAD
  NeuroX.05[[i]] <- br$NeuroX.FDR..05
  NeuroX.GWS[[i]] <- br$NeuroX.GWS
  ALSOD[[i]] <- br$ALSOD
  Cirulli[[i]] <- br$Carulli
}


# ######## EXAC ##########
# #Run for loop to extract each row into a separate list
# Exac_pval <- list()
# 
# #Collect enrichments for each list 
# for (i in 1:length(Exac)){
#   df <- as.data.frame(Exac[i
#                            ])
#   df <- as.data.frame(t(df))
#   row.names(df) <- i
#   
#   Exac_pval[[i]] <- df$BHadjP.value
# }

Exac_pval <- lapply(Exac,function(x)as.numeric(x[[3]]))
Exac_pval.df <- as.data.frame(Exac_pval)
Exac_pval.df <- t(Exac_pval.df)
rownames(Exac_pval.df) <- (1:nrow(Exac_pval.df))+2913

######## ALSOD ##########
ALSOD_pval <- lapply(ALSOD,function(x)as.numeric(x[[3]]))
ALSOD_pval.df <- as.data.frame(ALSOD_pval)
ALSOD_pval.df <- t(ALSOD_pval.df)
rownames(ALSOD_pval.df) <- (1:nrow(ALSOD_pval.df))+2913

######## Cirulli ##########
Cirulli_pval <- lapply(Cirulli,function(x)as.numeric(x[[3]]))
Cirulli_pval.df <- as.data.frame(Cirulli_pval)
Cirulli_pval.df <- t(Cirulli_pval.df)
rownames(Cirulli_pval.df) <- (1:nrow(Cirulli_pval.df))+2913

######## GeneCards AD ##########
GeneCards.AD_pval <- lapply(GeneCards.AD,function(x)as.numeric(x[[3]]))
GeneCards.AD_pval.df <- as.data.frame(GeneCards.AD_pval)
GeneCards.AD_pval.df <- t(GeneCards.AD_pval.df)
rownames(GeneCards.AD_pval.df) <- (1:nrow(GeneCards.AD_pval.df))+2913

######## GeneCards ALS ##########
GeneCards.ALS_pval <- lapply(GeneCards.ALS,function(x)as.numeric(x[[3]]))
GeneCards.ALS_pval.df <- as.data.frame(GeneCards.ALS_pval)
GeneCards.ALS_pval.df <- t(GeneCards.ALS_pval.df)
rownames(GeneCards.ALS_pval.df) <- (1:nrow(GeneCards.ALS_pval.df))+2913

######## GWASCentral ALS ##########
GWASCentral.ALS_pval <- lapply(GWASCentral.ALS,function(x)as.numeric(x[[3]]))
GWASCentral.ALS_pval.df <- as.data.frame(GWASCentral.ALS_pval)
GWASCentral.ALS_pval.df <- t(GWASCentral.ALS_pval.df)
rownames(GWASCentral.ALS_pval.df) <- (1:nrow(GWASCentral.ALS_pval.df))+2913

######## GWASCentral AD ##########
GWASCentral.AD_pval <- lapply(GWASCentral.AD,function(x)as.numeric(x[[3]]))
GWASCentral.AD_pval.df <- as.data.frame(GWASCentral.AD_pval)
GWASCentral.AD_pval.df <- t(GWASCentral.AD_pval.df)
rownames(GWASCentral.AD_pval.df) <- (1:nrow(GWASCentral.AD_pval.df))+2913

######## NeuroX.05 ##########
NeuroX.05_pval <- lapply(NeuroX.05,function(x)as.numeric(x[[3]]))
NeuroX.05_pval.df <- as.data.frame(NeuroX.05_pval)
NeuroX.05_pval.df <- t(NeuroX.05_pval.df)
rownames(NeuroX.05_pval.df) <- (1:nrow(NeuroX.05_pval.df))+2913

######## NeuroX.GWS ##########
NeuroX.GWS_pval <- lapply(NeuroX.GWS,function(x)as.numeric(x[[3]]))
NeuroX.GWS_pval.df <- as.data.frame(NeuroX.GWS_pval)
NeuroX.GWS_pval.df <- t(NeuroX.GWS_pval.df)
rownames(NeuroX.GWS_pval.df) <- (1:nrow(NeuroX.GWS_pval.df))+2913

######## Pasterkamp ##########
Pasterkamp_pval <- lapply(Pasterkamp,function(x)as.numeric(x[[3]]))
Pasterkamp_pval.df <- as.data.frame(Pasterkamp_pval)
Pasterkamp_pval.df <- t(Pasterkamp_pval.df)
rownames(Pasterkamp_pval.df) <- (1:nrow(Pasterkamp_pval.df))+2913

######## Taylor ##########
Taylor_pval <- lapply(Taylor,function(x)as.numeric(x[[3]]))
Taylor_pval.df <- as.data.frame(Taylor_pval)
Taylor_pval.df <- t(Taylor_pval.df)
rownames(Taylor_pval.df) <- (1:nrow(Taylor_pval.df))+2913



####################################################################################################################################################
####################################################################################################################################################
library(hgu133plus2.db)

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

for (i in 1:8000){
  C9_int <- c9_gene[1:i]
  CH_int <- ch_gene[1:i]
  sals_int <- sals_gene[1:i]
  ftld_int <- ftld_gene[1:i]
  vcp_int <- vcp_gene[1:i]
  pet_int <- pet_gene[1:i]
  rav_int <- rav_gene[1:i]
  List[[i]] <- Reduce(intersect, list(C9_int, CH_int, sals_int, ftld_int, vcp_int, pet_int, rav_int))
}

#Load file with all genes
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]
sym.genes <- t(sym.genes)
allgenes <- sym.genes[!duplicated(sym.genes),]

#Remove list elements with less than 5 genes (to aid calculations)
List_5 <- List[lengths(List) > 4]
#Leaves final 5087 elements (elements 1:2913 removed)

#Create new empty list
enrich_result <- list()

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/")
S <- read.table(file = "TDP-43_PPIgenes.txt")
s <- S$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/")
o <- read.table(file = "OneBenchmarkList.txt")
o <- o$V1




for (i in 1:length(List_5)){
  write.table(List_5[i], "benchmark_genelist.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  bgene <- read.table("benchmark_genelist.txt")
  bgene <- bgene$V1
  
  ur.list <- bgene
  int.list <- s

  #How many test geneset genes contain snps
  x.in <- length (which(ur.list %in% int.list)) 

  #how many do not
  x.out <- length(ur.list) - x.in

  #total number of snp genes
  tot.in <- length(int.list)

  #total number of all genes
  tot.out <- length(allgenes)-length(tot.in)


  #create count matrix
  counts <- matrix (nrow=2, ncol=2)
  counts [1,] <- c(x.in, tot.in)
  counts [2,] <- c(x.out, tot.out)


  #Conduct fisher's exact test for count data
  a5 <-fisher.test (counts)

  enrich_result[[i]] <- a5$p.value
}

enrich_result_df <- t(as.data.frame(enrich_result))
rownames(enrich_result_df) <- (1:nrow(enrich_result_df))+2913
ERpadj <- p.adjust(enrich_result_df, method = "BH")
ERpadj <- as.data.frame(ERpadj)
rownames(ERpadj) <- (1:nrow(ERpadj))+2913

cat(List[[6732]], sep = "\n")
length(List[[6732]])

write.table(List[[6732]], "LIST_6732.txt", row.names = F, col.names = F, quote = F)

logdf<- log(ERpadj)
logdf<- as.data.frame(logdf)
logdf[,2] <- rownames(logdf)
logdf[,3] <- logdf$ERpadj
logdf[,1] <- NULL
plot(logdf, xlab = "Threshold", ylab = "log P-Value", main = "Log P Value with increasing threshold")

int <- intersect(x = List[[6732]], y = s)
