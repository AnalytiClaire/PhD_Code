
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression")
C9 <- read.csv("C9rankeduniqueresult.csv")
sals <- read.csv("sALSrankeduniqueresult.csv")
ftld <- read.csv("FTLDrankeduniqueresult.csv")
vcp <- read.csv("VCPrankeduniqueresult.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_results_keepfiltering.csv")
pet <- pet[!duplicated(pet$hgnc_symbol),]
rav <- read.csv("RAV_results_keepfiltering.csv")
rav <- rav[!duplicated(rav$hgnc_symbol),]



setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NetworkAnalyst")
C9nug <- read.csv("C9_nugget_WC.csv", row.names = 1)
C9nug_p <- C9nug[,4:11]
C9nug_c <- C9nug[,1:3]
salsnug <- read.csv("sals_nugget_WC.csv", row.names = 1)
salsnug_p <- salsnug[,4:10]
salsnug_c <- salsnug[,1:3]
ftldnug <- read.csv("ftld_nugget_WC.csv", row.names = 1)
ftldnug_p <- ftldnug[,9:24]
ftldnug_c <- ftldnug[,1:8]
vcpnug <- read.csv("vcp_nugget_WC.csv", row.names = 1)
vcpnug_p <- vcpnug[,4:10]
vcpnug_c <- vcpnug[,1:3]
petnug <- read.csv("pet_nugget_WC.csv", row.names = 1)
petnug_p <- petnug[,10:26]
petnug_c <- petnug[,1:9]
ravnug <- read.csv("rav_nugget_WC.csv", row.names = 1)
ravnug_p <- ravnug[,9:21]
ravnug_c <- ravnug[,1:8]


genelist <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/TFRegression/GeneXplainTFs.txt")

subsetC9 <- subset(C9, C9$Gene.Symbol %in% genelist, drop = TRUE)
subsetC9 <- subsetC9[!duplicated(subsetC9$Gene.Symbol),]
rownames(subsetC9) <- subsetC9$Gene.Symbol
subsetC9 <- subsetC9[order(row.names(subsetC9)),]

subsetsals <- subset(sals, sals$Gene.Symbol %in% genelist, drop = TRUE)
subsetsals <- subsetsals[!duplicated(subsetsals$Gene.Symbol),]
rownames(subsetsals) <- subsetsals$Gene.Symbol
subsetsals <- subsetsals[order(row.names(subsetsals)),]

subsetftld <- subset(ftld, ftld$Gene.Symbol %in% genelist, drop = TRUE)
subsetftld <- subsetftld[!duplicated(subsetftld$Gene.Symbol),]
rownames(subsetftld) <- subsetftld$Gene.Symbol
subsetftld <- subsetftld[order(row.names(subsetftld)),]

subsetvcp <- subset(vcp, vcp$Gene.Symbol %in% genelist, drop = TRUE)
subsetvcp <- subsetvcp[!duplicated(subsetvcp$Gene.Symbol),]
rownames(subsetvcp) <- subsetvcp$Gene.Symbol
subsetvcp <- subsetvcp[order(row.names(subsetvcp)),]

subsetpet <- subset(pet, pet$hgnc_symbol %in% genelist, drop = TRUE)
subsetpet <- subsetpet[!duplicated(subsetpet$hgnc_symbol),]
rownames(subsetpet) <- subsetpet$Gene.Symbol
subsetpet <- subsetpet[order(subsetpet$hgnc_symbol),]
rownames(subsetpet) <- subsetpet$hgnc_symbol

subsetrav <- subset(rav, rav$hgnc_symbol %in% genelist, drop = TRUE)
subsetrav <- subsetrav[!duplicated(subsetrav$hgnc_symbol),]
rownames(subsetrav) <- subsetrav$Gene.Symbol
subsetrav <- subsetrav[order(subsetrav$hgnc_symbol),]
rownames(subsetrav) <- subsetrav$hgnc_symbol


C9_p <- subsetC9[,24:31]
sals_p <- subsetsals[,24:30]
ftld_p <- subsetftld[,28:43]
sfltd_p <- subsetftld[,34:43]
GRNfltd_p <- subsetftld[,28:33]
vcp_p <- subsetvcp[,24:30]
pet_p <- subsetpet[19:35]
petC9_p <- subsetpet[19:25]
petsals_p <- subsetpet[26:35]
rav_p <- subsetrav[,18:30]


tC9_p <- t(C9_p)
tC9nug <- t(C9nug_p)

plot(tC9_p[,1], tC9nug[,1])
x <- vector()
z <- vector()
result <- matrix(0,72,2)

for (i in 1:ncol(tC9_p)){
  for (j in 1:ncol(tC9nug)){
    plot(tC9_p[,i], tC9nug[,j], main = colnames(tC9nug[,i]))
    y <- cor.test(tC9_p[,i], tC9nug[,j])
    x[[i]]<- y$statistic
    z[[i]] <- y$p.value
  }
  result[,i] <- x
  result[,i+1] <- z
}



cor.test(tC9_p[,i], tC9nug[,i])

