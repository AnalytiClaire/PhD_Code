# # Example script for hyperPathway
# 
# library(pathprint)
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GSEA/PCxN Example/probesets/")
# 
# #Load individual gene names for each significance threshold
# A <- read.table(file = "threegenes.txt")
# a <- A$V1
# 
# B <- read.table(file = "fourgenes.txt")
# b <- B$V1
# 
# C <- read.table(file = "fivegenes.txt")
# c <- C$V1
# 
# D <- read.table(file = "sixgenes")
# d <- D$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS")
# 
# E <- read.table(file = "signif.snp.NeuroX.txt")
# e <- E$V1
# 
# F <- read.table(file = "signif.snp.NeuroX.p5E08.txt")
# f <- F$V1
# 
# setwd(dir = "/Users/clairegreen/Desktop/")
# 
# G <- read.table(file = "TDP-43DEGs.txt")
# g <- G$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS")
# H <- read.table(file = "signif.snp.AD.GWASCentralp5E08.txt")
# h <- H$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
# K <- read.table(file = "subnet.28.GM.txt")
# k <- K$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/ExAC/")
# L <- read.table(file = "exac.pli.0.95.txt")
# l <- L$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
# M <- read.table(file = "Cirulli.txt")
# m <- M$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
# N <- read.table(file = "GeneCardsAD.txt")
# n <- N$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
# O <- read.table(file = "GeneCardsALS.txt")
# o <- O$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/QQ/Test3/")
# P <- read.table(file = "genemania-genes.txt")
# p <- P$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
# Q <- read.table(file = "Pasterkamp_TDP43.txt")
# q <- Q$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
# R <- read.table(file = "Taylor_TDP43.txt")
# r <- R$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")
# S <- read.table(file = "C9_genenames.txt")
# s <- S$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/6500")
# U <- read.table("6500 + 200.txt")
# u <- U$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression")
# V <- read.table("OneBenchmarkList.txt")
# V <- V$V1
# 
# setwd(dir = "/Users/clairegreen/Desktop/")
# List5996 <- read.table("List5996.txt")
# List5996 <- List5996$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression")
# List6732 <- read.table("LIST_6732.txt")
# List6732 <- List6732$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
# up <- read.table("intersect_up_1.txt")
# up <- up$V1
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
# down <- read.table("intersect_down_1.txt")
# down <- down$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
up_filter <- read.table("Filtered_up_genes.txt")
up_filter <- up_filter$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
down_filter <-read.table("Filtered_down_genes.txt")
down_filter <- down_filter$V1

updown <- c(as.character(up_filter), as.character(down_filter))

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")
degppi <- read.table("DEG_PPI_Genes_nofib.txt")
degppi <- degppi$V1

Disgenes <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")

nugget <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NuggetGenes.txt")

wholenetwork_nug <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/Whole Network default node.csv")

# setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint 25.04.16/hyperPathway/All.Pathways/PathprintPathways(29)/")

# Z <- read.csv(file = "pathprintgenes.csv", na.strings = c("", "NA)"))
# Z <- as.list(Z)
# Z <- lapply(Z, function(x) x[!is.na(x)])

# z <- read.csv(file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/PathprintThreshold/GeneList_200.csv", na.strings = c("", "NA)"))
# z <- as.list(z)
# z <- lapply(z, function(x) x[!is.na(x)])

# SNF <- read.csv(file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/SimilarityNetworkFusion/microarray/Colourlist.csv", na.strings = c("", "NA)"))
# SNF <- as.list(SNF)
# SNF <- lapply(SNF, function(x) x[!is.na(x)])


# #Pathprint gene lists
# load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint_biomart.RData")

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

#Read in geneset
# setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
# Y <- read.csv(file = "AllgenesNO.csv", na.strings = c("", "NA)"))
# Y <- as.list(Y)
# Y<- lapply(Y, function(x) x[!is.na(x)])

# setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/LF_miRNA/")
# X <- read.csv(file = "miRNAtargetGenes.csv", na.strings = c("", "NA)"))
# X <- as.list(X)
# X<- lapply(X, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/")
A <- read.csv(file = "BenchmarkGenes2.csv", na.strings = c("", "NA)"))
A <- as.list(A)
A <- lapply(A, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/")
B <- read.csv(file = "Benchmarks3.csv", na.strings = c("", "NA)"))
B <- as.list(B)
B <- lapply(B, function(x) x[!is.na(x)])

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







# run script
pathwayEnrichment <- hyperPathway(
                genelist = deg,
								geneset = B,
								Nchip = length(allgenes)
							 )
setwd (dir = "/Users/clairegreen/Documents/PhD/MyPapers/TDP-43_1/")
write.csv(pathwayEnrichment, file = "BenchmarkEnrich_upDEG.csv")

nug <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/TDP_fixed_NuggetGenes.txt")
deg <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_upanddown.txt")

cat(intersect(nug, deg), sep = ", ")















o2 <- Reduce(intersect, list(up_filter, B$Disgenes))
cat(o2, sep = ", ")
write.csv(o2, file = "intersect.csv")

o2 <- Reduce(intersect, list(up_filter, B$Disgenes))
cat(o2, sep = ", ")
write.csv(o2, file = "intersect.csv")

o2 <- Reduce(intersect, list(up_filter, B$Disgenes))
cat(o2, sep = ", ")
write.csv(o2, file = "intersect.csv")



#As for loop
df <- data.frame(matrix(ncol=length(B), nrow=length(SNF)))
for (i in 1:length(B)){
  pathwayEnrichment <- hyperPathway(
    genelist = B[[i]],
    geneset = SNF,
    Nchip = length(allgenes)
  )
  df[,i] <- pathwayEnrichment$`BHadjP-value`
  colnames(df)[i] <- names(B[i])
  rownames(df) <- names(SNF)
}

write.csv(t(df), "BenchmarkEnrichmentTable(t).csv")
t(df)
