#Set working directory and load samples as lists where there are no NAs or empty values
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/")

C9Mod <- read.csv(file = "C9WGCNA.csv", na.strings = c("", "NA)"))
C9Mod <- as.list(C9Mod)
C9Mod<- lapply(C9Mod, function(x) x[!is.na(x)])

CHMod <- read.csv(file = "CHWGCNA.csv", na.strings = c("", "NA)"))
CHMod <- as.list(CHMod)
CHMod<- lapply(CHMod, function(x) x[!is.na(x)])

sALSMod <- read.csv(file = "sALSWGCNA.csv", na.strings = c("", "NA)"))
sALSMod <- as.list(sALSMod)
sALSMod<- lapply(sALSMod, function(x) x[!is.na(x)])

FTLDMod <- read.csv(file = "FTLDWGCNA.csv", na.strings = c("", "NA)"))
FTLDMod <- as.list(FTLDMod)
FTLDMod<- lapply(FTLDMod, function(x) x[!is.na(x)])

VCPMod <- read.csv(file = "VCPWGCNA.csv", na.strings = c("", "NA)"))
VCPMod <- as.list(VCPMod)
VCPMod<- lapply(VCPMod, function(x) x[!is.na(x)])

#Select two datasets to compare
geneset1 <- C9Mod
geneset2 <- sALSMod

#Create empty data frames to hold the output
intersect <- as.data.frame(matrix(nrow = length(geneset1), ncol = length(geneset2)))
rownames(intersect) <- names(geneset1)
colnames(intersect) <- names(geneset2)

perc_1 <- as.data.frame(matrix(nrow = length(geneset1), ncol = length(geneset2)))
rownames(perc_1) <- names(geneset1)
colnames(perc_1) <- names(geneset2)

perc_2 <- as.data.frame(matrix(nrow = length(geneset1), ncol = length(geneset2)))
rownames(perc_2) <- names(geneset1)
colnames(perc_2) <- names(geneset2)

#Do cross comparison of all list members from the two lists
for (i in 1:length(geneset1)) {
  for (j in 1:length(geneset2)){
  int <- length(intersect(geneset1[[i]], geneset2[[j]]))
  intersect[i,j] <- int
  
  perc1 <- (int/length(geneset1[[i]]))*100
  perc_1[i,j] <- perc1
  
  perc2 <- (int/length(geneset2[[j]]))*100
  perc_2[i,j] <- perc2
  
  }
}

#Save in directory
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/WGCNA_overlap/")
write.csv(intersect, file = "C9_sals_intersect.csv")
write.csv(perc_1, file = "C9_sals_gs1perc.csv")
write.csv(perc_2, file = "C9_sals_gs2perc.csv")