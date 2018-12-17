# set working directory
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HoffmanMuscle/")

#Load individual gene names for each significance threshold
A <- read.table(file = "MD Genes.txt")
a <- A$V1

B <- read.table(file = "CongenitalMDGenes.txt") ###
b <- B$V1

C <- read.table(file = "CongenitalMyopathies.txt")
c <- C$V1

D <- read.table(file = "OtherMyopathies.txt")
d <- D$V1

E <- read.table(file = "MotorNeuronDiseases.txt") ###
e <- E$V1

F <- read.table(file = "DistalMyopathies.txt")
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

x <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HoffmanMuscle/SigGenes.01.txt")
ur.list <- x$V1
int.list <- f

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