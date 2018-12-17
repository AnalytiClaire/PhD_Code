## RPT for common up and common down ###

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
DIJ <- read.csv("DIJfilteredresult.csv")
DUM <- read.csv("DUM_UniqueGene_DESeq2.csv")
FFR <- read.csv("FFRfilteredresult.csv")
LEW <- read.csv("LEWfilteredresult.csv")
MID1 <- read.csv("MID1filteredresult.csv")
MID2 <- read.csv("MID2filteredresult.csv")
MID3 <- read.csv("MID3filteredresult.csv")
MID4 <- read.csv("MID4filteredresult.csv")
MOR.FC <- read.csv("MOR.FCfilteredresult.csv")
MOR.SN <- read.csv("MOR.SNfilteredresult.csv")
BOT <- read.csv("BOTrankeduniqueresult.csv")
BOT2 <- read.csv("BOT2rankeduniqueresult.csv")


DIJgene <- as.character(DIJ$Gene.Symbol)
DUMgene <- as.character(DUM$hgnc_symbol)
FFRgene <- as.character(FFR$Gene.Symbol)
LEWgene <- as.character(LEW$Gene.Symbol)
MID1gene <- as.character(MID1$Gene.Symbol)
MID2gene <- as.character(MID2$Gene.Symbol)
MID3gene <- as.character(MID3$Gene.Symbol)
MID4gene <- as.character(MID4$Gene.Symbol)
MOR.FCgene <- as.character(MOR.FC$Gene.Symbol)
MOR.SNgene <- as.character(MOR.SN$Gene.Symbol)
BOTgene <- as.character(BOT$Gene.Symbol)
BOT2gene <- as.character(BOT2$Gene.Symbol)





m = 100000
r <- matrix(0, m, 3)

for (i in 1:m){

  DIJrand <- sample(DIJgene)
  upDIJ <- DIJrand[1:9728]
  downDIJ <- DIJrand[9729:23346]
  
  DUMrand <- sample(DUMgene)
  upDUM <- DUMrand[1:7854]
  downDUM <- DUMrand[7855:20490]
  
  FFRrand <- sample(FFRgene)
  upFFR <- FFRrand[1:6582]
  downFFR <- FFRrand[6583:13433]
  
  LEWrand <- sample(LEWgene)
  upLEW <- LEWrand[1:7664]
  downLEW <- LEWrand[7665:13433]
  
  MID1rand <- sample(MID1gene)
  upMID1 <- MID1rand[1:6489]
  downMID1 <- MID1rand[6490:13433]
  
  MID2rand <- sample(MID2gene)
  upMID2 <- MID2rand[1:7529]
  downMID2 <- MID2rand[7530:13433]
  
  MID3rand <- sample(MID3gene)
  upMID3 <- MID3rand[1:8331]
  downMID3 <- MID3rand[8332:13433]
  
  MID4rand <- sample(MID4gene)
  upMID4 <- MID4rand[1:7174]
  downMID4 <- MID4rand[7175:13433]
  
  MOR.FCrand <- sample(MOR.FCgene)
  upMOR.FC <- MOR.FCrand[1:7829]
  downMOR.FC <- MOR.FCrand[7830:13433]
  
  MOR.SNrand <- sample(MOR.SNgene)
  upMOR.SN <- MOR.SNrand[1:7944]
  downMOR.SN <- MOR.SNrand[7945:13433]
  
  BOTrand <- sample(BOTgene)
  upBOT <- BOTrand[1:15045]
  downBOT <- BOTrand[15046:28262]
  
  BOT2rand <- sample(BOT2gene)
  upBOT2 <- BOT2rand[1:15592]
  downBOT2 <- BOT2rand[15593:28262]
  

  
  INTUP <- Reduce(intersect, list(upDIJ, upDUM, upFFR, upLEW, upMID1, upMID2, 
                                  upMID3, upMID4, upMOR.FC, upMOR.SN, upBOT, upBOT2))
  r[i,1] <- length(INTUP)
  

  INTDOWN <- Reduce(intersect, list(downDIJ, downDUM, downFFR, downLEW, downMID1, downMID2, 
                                  downMID3, downMID4, downMOR.FC, downMOR.SN, downBOT, downBOT2))
  r[i,2] <- length(INTDOWN)
  r[i,3] <- sum(length(INTUP) + length(INTDOWN)) 
}

setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
r <- read.csv("UpDownRPT.csv")

r <- as.data.frame(r)
expup <- 175
expdown <- 189
exptotal <- 364

# write.csv(r, "UpDownRPT.csv")

testup <- which(r$V1 >= expup) 
resultup <- sum((length(testup)+1))/(m+1) # calculate P value
resultup
mean <- mean(r$V1)
mean
range <- range(r$V1)
range

hist(r$V1, 
     xlim = range(0:200), 
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
     xlim = range(0:200), 
     main = NULL, 
     xlab = "Number of Common Downregulated DEGs")
abline(v = 189, col = "red", lwd = 2)

testtotal <- which(r$V3 >= exptotal) 
resulttotal <- sum((length(testtotal)+1))/(m+1) # calculate P value
resulttotal
mean <- mean(r$V3)
mean
range <- range(r$V3)
range

hist(r$V3, 
     xlim = range(0:260), 
     main = NULL, 
     xlab = "Number of Common DEGs")
abline(v = 243, col = "red", lwd = 2)



table <- data.frame(NumOverTest = length(test1),
                    Pval = result,
                    mean = mean,
                    range = range)

