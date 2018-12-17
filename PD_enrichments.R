#Fishers Exact Test


setwd
PARA_DEGs <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PARA/GeneExpression/Fixed_JENRA/allDEGs_OA_Blood_remove.txt")
PARA_PPI <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PARA/GeneExpression/Fixed_JENRA/DEG_PPI_Genes.txt")
PARA_Nug <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PARA/GeneExpression/Fixed_JENRA/JENRA_Nug_Genes.txt")
RAGenes <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/RAGenes.txt")
OAGenes <- readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/OAGenes.txt")
celltype <- read.csv("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv")
DEG <- readLines("../GeneExpression/allDEGs_SJO_SLE_ALO.txt")
Okada <- readLines("../../Okada_MetaGWAS_RA.txt")

x <- read.csv("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/ComCoexpr/RAOA_0point3_NoWAL.OAPAD.node.csv")



TDP43PPI <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes_nofib.txt")
TDP43gene <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")
TDPmodule <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/TDP_fixed_NuggetGenes.txt")


## PD 
PD_bench <- read.csv("/Users/clairegreen/Documents/PhD/Parkinsons/PD_Benchmarks.csv")
PD_gwas <- as.character(PD_bench$PD.GWAS)
PD_lrrk2ppi <- as.character(PD_bench$LRRK2.PPI)
PD_DEGs <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/ALS_sfblood_ALLgenes.txt")
PD_PPI <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PPIGenes.txt")
PD_Nug <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PD_Nuggenes_confambloodremoved.txt")

#### RUN FISHER'S EXACT TEST ###


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

# 
# ur.list <- PD_DEGs
# int.list <- PD_gwas
# 

# #How many test geneset genes contain snps
# x.in <- length (which(ur.list %in% int.list)) 
# #how many do not
# x.out <- length(ur.list) - x.in
# #total number of snp genes
# tot.in <- length(int.list)
# #total number of all genes
# tot.out <- length(allgenes)-length(tot.in)
# 
# 
# #create count matrix
# counts <- matrix (nrow=2, ncol=2)
# counts [1,] <- c(x.in, tot.in)
# counts [2,] <- c(x.out, tot.out)
# 
# #Conduct fisher's exact test for count data
# a5 <-fisher.test(counts)
# enrich <- a5$p
# print(enrich)




######

#### See enrichment with TDP-43 benchmark genes ###

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

#####
hyperPathway <-
  function(genelist, geneset, Nchip) # produces a list of pathways enrichments calculated using the hypergeometric distribution
  {
    # Calculate p-values

    Nsig <- length(genelist)
    hyper<-as.data.frame(matrix(nrow = length(geneset), ncol = 1))
    rownames(hyper) <- names(geneset)
    colnames(hyper) <- c("p-value")
    # determine p-value using the hypergeometric distribution, setting p = 1 if overlap = 0

    for (i in 1:length(geneset)){
      if (length(intersect(genelist, unlist(geneset[i]))) < 1){
        hyper[i,1]<-1}

      else if (length(intersect(genelist, unlist(geneset[i]))) > 0){
        hyper[i,1]<-phyper(length(intersect(genelist, unlist(geneset[i]))), Nsig, Nchip-Nsig, length(unlist(geneset[i])), lower.tail = FALSE)
      }
    }

    # adjust for multiple testing using Benjamini & Hochberg (fdr) correction
    hyper[,2]<-p.adjust(hyper[,1], method = "BH")
    # Obtain list genes for each pathway

    overlap <- vector("list", 0)
    for (i in 1:length(geneset)){
      temp.overlap <- list(intersect(genelist, unlist(geneset[[i]])))
      overlap <- append(overlap, temp.overlap)			}

    names(overlap)<-rownames(hyper)

    # Count number of list genes and total genes in each pathway

    for (i in 1:length(geneset)){
      hyper[i,3] <- length(overlap[[i]])
      hyper[i,4] <- length(geneset[[i]])
    }

    hyper[,5]<-rownames(hyper)
    hyper<-cbind((1:length(hyper[,1])), hyper)
    colnames(hyper)<-c("ID", "P-value", "BHadjP-value", "nGenes", "nPathway", "Name")
    return(hyper)
  } # end of hyperPathway
#####

setwd (dir = "/Users/clairegreen/Documents/PhD/Parkinsons/")
B <- read.csv(file = "PD_Benchmarks.csv", na.strings = c("", "NA)"))
B <- as.list(B)
B <- lapply(B, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/Arthritis/")
C <- read.csv(file = "Arthritis_Benchmarks.csv", na.strings = c("", "NA)"))
C <- as.list(C)
C <- lapply(C, function(x) x[!is.na(x)])
C <- C[1:2]


setwd (dir = "/Users/clairegreen/Documents/PhD/Arthritis/GWASCatalog/")
D <- read.csv(file = "GWAS_Catalog_ReportedGene.csv", na.strings = c("", "NA)"))
D <- as.list(D)
D <- lapply(D, function(x) x[!is.na(x)])


setwd (dir = "/Users/clairegreen/Documents/PhD/Arthritis/GWASCatalog/")
E <- read.csv(file = "GWAS_Catalog_ReportedGene_GWS.csv", na.strings = c("", "NA)"))
E <- as.list(E)
E <- lapply(E, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF")
F <- read.csv(file = "TDP_GWAS.csv", na.strings = c("", "NA)"))
F <- as.list(F)
F <- lapply(F, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results")
G <- read.csv(file = "Benchmarks3.csv", na.strings = c("", "NA)"))
G <- as.list(G)
G <- lapply(G, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results")
H <- read.csv(file = "BenchmarkGenes2.csv", na.strings = c("", "NA)"))
H <- as.list(H)
H <- lapply(H, function(x) x[!is.na(x)])


# run script
pathwayEnrichment <- hyperPathway(
  genelist = PARA_Nug,
  geneset = C,
  Nchip = length(allgenes)
)

cat(intersect(F$AD_GWAS, TDPmodule), sep = ", ")

cat(PD_DEGs, sep = "\n")
