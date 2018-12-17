# set working directory
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GSEA/PCxN Example/probesets/")

#Load individual gene names for each significance threshold
A <- read.table(file = "threegenes.txt")
a <- A$V1

B <- read.table(file = "fourgenes.txt")
b <- B$V1

C <- read.table(file = "fivegenes.txt")
c <- C$V1

D <- read.table(file = "sixgenes")
d <- D$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS")

E <- read.table(file = "signif.snp.NeuroX.txt")
e <- E$V1

F <- read.table(file = "signif.snp.NeuroX.p5E08.txt")
f <- F$V1

setwd(dir = "/Users/clairegreen/Desktop/")

G <- read.table(file = "Claire_intersect.txt")
g <- G$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS")
H <- read.table(file = "signif.snp.AD.GWASCentralp5E08.txt")
h <- H$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
K <- read.table(file = "subnet.28.GM.txt")
k <- K$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/ExAC/")
L <- read.table(file = "exac.pli.0.95.txt")
l <- L$V1


setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
M <- read.table(file = "Carulli_List.txt")
m <- M$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
N <- read.table(file = "GeneCardsAD.txt")
n <- N$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
O <- read.table(file = "GeneCardsALS.txt")
o <- O$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/QQ/Test3/")
P <- read.table(file = "genemania-genes.txt")
p <- P$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
Q <- read.table(file = "Pasterkamp_TDP43.txt")
q <- Q$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
R <- read.table(file = "Taylor_TDP43.txt")
r <- R$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/")
S <- read.table(file = "OneBenchmarkList.txt")
s <- S$V1


setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint 25.04.16/FishersExact/FE.All.Pathways/FE.PathprintPathways(29)/")
x <- read.table(file = "FE.pathprintgenes.csv")
x <- x$V1


setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/")
LC_TRL <- read.table(file = "LC/Q331K GRASPS_BitSeq_EdgeR_DE transcripts.txt")
LC_TRL <- LC_TRL$V3



Z <- read.csv(file = "FE.pathprintgenes.csv", na.strings = c("", "NA)"))
Z <- as.list(Z)


#Intersect
overlap <- Reduce(intersect, list(o2, o3))
print(overlap)



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


#calculate counts

# n <- max(length(a), length(b), length(c), length(d))
# length(a) <- n
# length(b) <- n
# length(c) <- n
# length(d) <- n
# 
# x <- cbind(a,b,c,d)
# 
# for (i in x[,1:4]) {

x <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/highly2014/highly5500.txt")
ur.list <- X
int.list <- Sig.padj.both

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
enrich <- a5$p
print(enrich)



o1 <- Reduce(intersect, list(ur.list, int.list))
print(o1)

o2 <- Reduce(intersect, list(f, Z$Signaling.by.Insulin.receptor..Reactome.))
print(o2)

