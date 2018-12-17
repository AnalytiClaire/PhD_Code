## RPT for common up and common down ###
C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
sals <- read.csv("sals_unique.csv")
sals <- sals[order(sals$P.Value),]
ftld <- read.csv("ftld_unique.csv")
ftld <- ftld[order(ftld$P.Value),]
vcp <- read.csv("vcp_unique.csv")
vcp <- vcp[order(vcp$P.Value),]
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

write.csv(r, "UpDownRPT.csv", row.names = F)