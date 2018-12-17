###Random permutations for Differential Expression##
options(scipen=999)

DEGPPI <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes_nofib.txt")
DEGs <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_upanddown.txt")
ALSOD <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/ALSoDgenes.txt")
Disgenes <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")
Taylor <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/Taylor_TDP43.txt")
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


setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/")
C9_cor.5 <- read.csv("C9_cor.5_cytoscape.csv")
sALS_cor.5 <-read.csv("sALS_cor.5_cytoscape.csv")
FTLD_cor.5 <-read.csv("FTLD_cor.5_cytoscape.csv")
VCP_cor.5 <-read.csv("VCP_cor.5_cytoscape.csv")
PET_cor.5 <-read.csv("PET_cor.5_cytoscape.csv")
RAV_cor.5 <-read.csv("RAV_cor.5_cytoscape.csv")

#Find length of overlap
length(intersect(nugget, Taylor))

#indicate the number of overlapping genes identified by DE analysis
test <- 4
samplenum <- 72
samplelist <- Taylor

m=100000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"


for (j in 1:m){
  PPIsample <- sample(allgenes, size=2550, replace=FALSE)
  modulesample <- sample(PPIsample, size = samplenum, replace = FALSE)
  random <- Reduce(intersect, list(modulesample, samplelist))
  r[j] <- length(random)
}

test1 <- which(r > test)  # count number of times r is larger than test value
result <- sum((length(test1)+1))/(m+1) # calculate P value
result
mean(r)
range(r)

###########################
#### label permutation ####
###########################

C9 <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/C9_cormatrix.csv", row.names = 1)
sALS <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/sALS_cormatrix.csv", row.names = 1)
FTLD <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/FTLD_cormatrix.csv", row.names = 1)
VCP <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/VCP_cormatrix.csv", row.names = 1)
PET <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PET_cormatrix.csv", row.names = 1)
RAV <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/RAV_cormatrix.csv", row.names = 1)

library(parallel)
library(BiocGenerics)
library(graph)
library(GraphAT)





C9_ordered <- C9[order(rownames(C9)),order(colnames(C9))]
sALS_ordered <- sALS[order(rownames(sALS)),order(colnames(sALS))]
FTLD_ordered <- FTLD[order(rownames(FTLD)),order(colnames(FTLD))]
VCP_ordered <- VCP[order(rownames(VCP)),order(colnames(VCP))]
PET_ordered <- PET[order(rownames(PET)),order(colnames(PET))]
RAV_ordered <- RAV[order(rownames(RAV)),order(colnames(RAV))]

C9_ordered <- as.matrix(C9_ordered)
sALS_ordered <- as.matrix(sALS_ordered)
FTLD_ordered <- as.matrix(FTLD_ordered)
VCP_ordered <- as.matrix(VCP_ordered)
PET_ordered <- as.matrix(PET_ordered)
RAV_ordered <- as.matrix(RAV_ordered)

x <- sample(rownames(C9_ordered))
key <- data.frame(Original = rownames(C9_ordered),
                  Random = x)

C9_permlabel <- C9_ordered
rownames(C9_permlabel) <- colnames(C9_permlabel) <- x
tri <- C9_permlabel
tri[lower.tri(tri, diag = TRUE)] <- NA
tri <- unmatrix(tri)
vec <- na.omit(tri)
vec <- as.data.frame(vec)
colnames(vec) <- "Rho"
cor.5.C9 <- subset(vec, vec$Rho > 0.5 | vec$Rho < -0.5)

sALS_permlabel <- sALS_ordered
rownames(sALS_permlabel) <- colnames(sALS_permlabel) <- x
tri <- sALS_permlabel
tri[lower.tri(tri, diag = TRUE)] <- NA
tri <- unmatrix(tri)
vec <- na.omit(tri)
vec <- as.data.frame(vec)
colnames(vec) <- "Rho"
cor.5.sALS <- subset(vec, vec$Rho > 0.5 | vec$Rho < -0.5)

FTLD_permlabel <- FTLD_ordered
rownames(FTLD_permlabel) <- colnames(FTLD_permlabel) <- x
tri <- FTLD_permlabel
tri[lower.tri(tri, diag = TRUE)] <- NA
tri <- unmatrix(tri)
vec <- na.omit(tri)
vec <- as.data.frame(vec)
colnames(vec) <- "Rho"
cor.5.FTLD <- subset(vec, vec$Rho > 0.5 | vec$Rho < -0.5)

VCP_permlabel <- VCP_ordered
rownames(VCP_permlabel) <- colnames(VCP_permlabel) <- x
tri <- VCP_permlabel
tri[lower.tri(tri, diag = TRUE)] <- NA
tri <- unmatrix(tri)
vec <- na.omit(tri)
vec <- as.data.frame(vec)
colnames(vec) <- "Rho"
cor.5.VCP <- subset(vec, vec$Rho > 0.5 | vec$Rho < -0.5)

PET_permlabel <- PET_ordered
rownames(PET_permlabel) <- colnames(PET_permlabel) <- x
tri <- PET_permlabel
tri[lower.tri(tri, diag = TRUE)] <- NA
tri <- unmatrix(tri)
vec <- na.omit(tri)
vec <- as.data.frame(vec)
colnames(vec) <- "Rho"
cor.5.PET <- subset(vec, vec$Rho > 0.5 | vec$Rho < -0.5)

RAV_permlabel <- RAV_ordered
rownames(RAV_permlabel) <- colnames(RAV_permlabel) <- x
tri <- RAV_permlabel
tri[lower.tri(tri, diag = TRUE)] <- NA
tri <- unmatrix(tri)
vec <- na.omit(tri)
vec <- as.data.frame(vec)
colnames(vec) <- "Rho"
cor.5.RAV <- subset(vec, vec$Rho > 0.5 | vec$Rho < -0.5)


Commonedge <- Reduce(intersect, list(rownames(cor.5.C9), rownames(cor.5.sALS), rownames(cor.5.FTLD),
                                  rownames(cor.5.VCP), rownames(cor.5.PET), rownames(cor.5.RAV)))


#Subset each dataset with these common names so they are all the same size
C9_CE <- subset(cor.5.C9, rownames(cor.5.C9) %in% Commonedge)
VCP_CE <- subset(cor.5.VCP, rownames(cor.5.VCP) %in% Commonedge)
FTLD_CE <- subset(cor.5.FTLD, rownames(cor.5.FTLD) %in% Commonedge)
sALS_CE <- subset(cor.5.sALS, rownames(cor.5.sALS) %in% Commonedge)
PET_CE <- subset(cor.5.PET, rownames(cor.5.PET) %in% Commonedge)
RAV_CE <- subset(cor.5.RAV, rownames(cor.5.RAV) %in% Commonedge)

m <- merge(C9_CE, sALS_CE, by = 0)
m <- merge(m, FTLD_CE, by.x = "Row.names", by.y = 0)
m <- merge(m, VCP_CE, by.x = "Row.names", by.y = 0)
m <- merge(m, PET_CE, by.x = "Row.names", by.y = 0)
m <- merge(m, RAV_CE, by.x = "Row.names", by.y = 0)

rownames(m) <- m[,1]
m[,1] <- NULL
colnames(m) <- c("C9", "sALS", "FTLD", "VCP", "PET", "RAV")
m_conserved <- m[apply(m, MARGIN = 1, function(x) all(x > 0)), ]





m <- Reduce(function(x,y) merge(x,y,0,all=T),list(C9_CE,VCP_CE,FTLD_CE,sALS_CE,PET_CE,RAV_CE))

vec$Genes <- rownames(vec)
vec$Gene1 <- as.character(lapply(strsplit(as.character(vec$Genes), "\\:"), "[", 2))
vec$Gene2 <- as.character(lapply(strsplit(as.character(vec$Genes), "\\:"), "[", 1))
