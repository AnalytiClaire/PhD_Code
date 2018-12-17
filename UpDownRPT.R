## RPT for common up and common down ###

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
sals <- read.csv("sals_unique.csv")
sals <- sals[order(sals$P.Value),]
ftld <- read.csv("ftld_unique.csv")
ftld <- ftld[order(ftld$P.Value),]
vcp <- read.csv("vcp_unique.csv")
vcp <- vcp[order(vcp$P.Value),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
rav <- read.csv("RAV_results_keepfiltering.csv")

m = 100000
r <- matrix(0, m, 3)

for (i in 1:m){
  #Sample from all genes "up" genes of the same size as experiment. This means the overlap is proportional.
  upC9 <- sample(C9$Gene.Symbol, size = 3788)
  upC9 <- as.vector(upC9)
  upSALS <- sample(sals$Gene.Symbol, size = 5905)
  upSALS <- as.vector(upSALS)
  upFTLD <- sample(ftld$Gene.Symbol, size = 4941)
  upFTLD <- as.vector(upFTLD)
  upVCP <- sample(vcp$Gene.Symbol, size = 8011)
  upVCP <- as.vector(upVCP)
  upPET <- sample(pet$hgnc_symbol, size = 9259)
  upPET <- as.vector(upPET)
  upRAV <- sample(rav$hgnc_symbol, size = 8028)
  upRAV <- as.vector(upRAV)
  
  INTUP <- Reduce(intersect, list(upC9, upSALS, upFTLD, upVCP, upPET, upRAV))
  r[i,1] <- length(INTUP)
  
  #### DOWN ####
  thresh <- -1
  
  downC9 <- subset(C9, !(C9$Gene.Symbol %in% upC9))
  downC9 <- downC9$Gene.Symbol
  downSALS <- subset(sals, !(sals$Gene.Symbol %in% upSALS))
  downSALS <- downSALS$Gene.Symbol
  downFTLD <- subset(ftld, !(ftld$Gene.Symbol %in% upFTLD))
  downFTLD <- downFTLD$Gene.Symbol
  downVCP <- subset(vcp, !(vcp$Gene.Symbol %in% upVCP))
  downVCP <- downVCP$Gene.Symbol
  downPET <- subset(pet, !(pet$hgnc_symbol %in% upPET))
  downPET <- downPET$hgnc_symbol
  downRAV <- subset(rav, !(rav$hgnc_symbol %in% upRAV))
  downRAV <- downRAV$hgnc_symbol
  
  INTDOWN <- Reduce(intersect, list(downC9, downSALS, downFTLD, downVCP, downPET, downRAV))
  r[i,2] <- length(INTDOWN)
  r[i,3] <- sum(length(INTUP) + length(INTDOWN)) 
}

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults")
r <- read.csv("UpDownRPT.csv")
r <- read.csv("UpDownRPT.csv")
expup <- 328
expdown <- 69
exptotal <- 397

testup <- which(r$V1 >= expup) 
resultup <- sum((length(testup)+1))/(m+1) # calculate P value
resultup
mean <- mean(r$V1)
mean
range <- range(r$V1)
range

hist(r$V1, 
     xlim = range(50:expup+30), 
     main = NULL, 
     xlab = "Number of Common Upregulated DEGs")
abline(v = expup, col = "red", lwd = 2)

testdown <- which(r$V2 >= expdown) 
resultdown <- sum((length(testdown)+1))/(m+1) # calculate P value
resultdown
mean <- mean(r$V2)
mean
range <- range(r$V2)
range

hist(r$V2, 
     xlim = range(0:80), 
     main = NULL, 
     xlab = "Number of Common Downregulated DEGs")
abline(v = expdown, col = "red", lwd = 2)

testtotal <- which(r$V3 >= exptotal) 
resulttotal <- sum((length(testtotal)+1))/(m+1) # calculate P value
resulttotal
mean <- mean(r$V3)
mean
range <- range(r$V3)
range

hist(r$V3, 
     xlim = range(80:exptotal+50), 
     main = NULL, 
     xlab = "Number of Common DEGs")
abline(v = exptotal, col = "red", lwd = 2)



table <- data.frame(NumOverTest = length(test1),
                    Pval = result,
                    mean = mean,
                    range = range)

