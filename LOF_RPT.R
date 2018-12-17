# GET ALL GENE NAMES
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

# FIND PROTEIN CODING GENES' START AND END POSITIONS
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot", "start_position", "end_position"), filters="hgnc_symbol", values=allgenes,  mart=mart)

# REMOVE NON-CODING GENES
mb <- subset(mart_back, mart_back$uniprotswissprot != "")
mb <- mb[!duplicated(mb[,1]),]

# FIND LENGTH OF GENE BY SUBTRATING START FROM END POSITION
mb$length <- mb$end_position - mb$start_position
hist(mb$length)

# READ IN EXPERIMENTAL GENE SET AND SUBSET
nug <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/TDP_NuggetGenes_genetics.txt")
mart_nug <- subset(mb, mb$hgnc_symbol %in% nug)

#MAKE SURE ALL GENES ARE PRESENT
setdiff(nug, mart_nug$hgnc_symbol)

# CONVERT LENGTHS TO A LIST
genelen <- mart_nug$length

# CREATE A LIST OF GENE LENGTHS WITHIN 1000BP EITHER SIDE OF TEST GENE LENGTHS
# SUBSET ALL GENES TO FIND GENES WITHIN THAT RANGE
# SAVE SYMBOLS TO A LIST
pop <- list()

for (i in 1:length(nug)){
  y <- seq(genelen[i]-1000, genelen[i]+1000, by=1)
  x <- subset(mb, mb$length %in% y)
  pop[[i]] <- x$hgnc_symbol
}

poplength <- vector()
for (i in 1:length(pop)){
  poplength[i] <- length(pop[[i]])
}

hist(poplength)

lof1 <- read.csv("All_data_ALS_1.csv")
lof2 <- read.csv("All_data_ALS_2.csv")
lof <- rbind(lof1, lof2)
lof_rare <- subset(lof, lof$is_lof == 1 & lof$aaf_exac_all < 0.01 & lof$aaf_1kg_all < 0.01)
#lof_rare <- subset(lof, lof$is_lof == 1 & lof$in_1kg == 0 & lof$in_exac == 0)

m = 1000
result_lof <- vector()
result_gene <- vector()

for (i in 1:m){
  
  selection <- vector()
    
  for (j in 1:length(pop)){
      selection[j] <- sample(pop[[j]])
    }
  
  rand_lof <- subset(lof_rare, lof_rare$gene %in% selection)
  result_lof[i] <- nrow(rand_lof)
  gene <- as.character(rand_lof$gene)
  gene <- unique(gene)
  result_gene[i] <- length(gene)
  
}

nug_lof <- read.csv("Nug_LOF_ProjectMinE.csv")

exp <- 77

test <- which(result_gene >= exp) 
res <- sum((length(test)+1))/(m+1) # calculate P value
res
mean <- mean(result)
mean
range <- range(result)
range

hist(result, 
     xlim = range(0:500), 
     main = NULL, 
     xlab = "Number of LoF Mutations")
abline(v = expup, col = "red", lwd = 2)