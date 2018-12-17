##### RNA seq analysis comparisons #####

RNA_deseq2 <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/Guillaume_HEK48h/DEseq2/CG_TRL_Gene_DESeq2.txt")
RNA_deseq2 <- RNA_deseq2$V1

RNA_all <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/Guillaume_HEK48h/Limma/TRL_Comparison.csv")
RNA_limma <- RNA_all$CG

# SA_TRL <- RNA_all$SA
# LC_TRL <- RNA_all$LC

RNA_EdgeR <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/Guillaume_HEK48h/EdgeR/TRL_diffexpr-genes.txt")
RNA_EdgeR <- RNA_EdgeR$V1

common_DEGs <- read.table("/Users/clairegreen/Desktop/TDP-43DEGs.txt")
common_DEGs <- common_DEGs$V1

overlap <- Reduce(intersect, list(LC_TRL, common_DEGs))
print(overlap)


library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)])
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

x <- LC_TRL
y <- common_DEGs


x.in <- length (which(x %in% y))
x.out <- length(x) - x.in
tot.in <- length (y)
tot.out <- length (sym.genes)

counts <- matrix (nrow=2, ncol=2)
counts [1,] <- c(x.in, tot.in)
counts [2,] <- c(x.out, tot.out)

a5 <-fisher.test (counts)
enrich <- a5$p
