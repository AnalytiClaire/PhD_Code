#Benchmarking for Paper

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

#Benchmarking Genes
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


genelist <- nugget

# run script
pathwayEnrichment <- hyperPathway(
  genelist = genelist,
  geneset = B,
  Nchip = length(allgenes)
)
setwd (dir = "/Users/clairegreen/Documents/PhD/MyPapers/TDP-43_1/")
write.csv(pathwayEnrichment, file = "BenchmarkEnrich_PPI.csv")




o1 <- Reduce(intersect, list(genelist, B$Disgenes))
cat(o1, sep = ", ")

o2 <- Reduce(intersect, list(genelist, B$Subnetwork28))
cat(o2, sep = ", ")

o3 <- Reduce(intersect, list(genelist, B$TDP.43.PPI))
cat(o3, sep = ", ")
