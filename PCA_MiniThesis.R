library(ggplot2)
library(gridExtra)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
C9 <- read.csv("C9rankeduniqueresult.csv")
C9pat <- C9[,52:59]
CH <- read.csv("CHrankeduniqueresult.csv")
CHpat <- CH[,55:57]
sALS <- read.csv("sALSrankeduniqueresult.csv")
sALSpat <- sALS[,52:58]
FTLD <- read.csv("FTLDrankeduniqueresult.csv")
FTLDpat <- FTLD[,57:72]
VCP <- read.csv("VCPrankeduniqueresult.csv")
VCPpat <- VCP[,52:58]
pet <- read.csv("PETNOrankeduniqueresult.csv")
Petpat <- pet[,52:59]
rav <- read.csv("RAVEXPRSrankeduniqueresult.csv")
Ravpat <- rav[,52:59]

#PCA on samples
pcaC9 <- prcomp(t(C9pat)) #run pca
pcaCH <- prcomp(t(CHpat)) #run pca
pcasALS <- prcomp(t(sALSpat)) #run pca
pcaFTLD <- prcomp(t(FTLDpat)) #run pca
pcaVCP <- prcomp(t(VCPpat)) #run pca

dev.off() #removes previous plot
plot(pcaC9$x[,1:2], pch=18, cex=1.25 , col="red",main = "PCA of Patient Microarray Datasets", 
     xlim=c(-20000, 20000), ylim=c(-20000, 20000))
points(pcaCH$x[,1:2], pch=18, cex=1.25 , col="orange")
points(pcasALS$x[,1:2], pch=18, cex=1.25 , col="green")
points(pcaFTLD$x[,1:2], pch=18, cex=1.25 , col="blue")
points(pcaVCP$x[,1:2], pch=18, cex=1.25 , col="magenta")

legend(16000, 20000 ,pch=15, legend=c("C9orf72", "CHMP2B", "sALS", "FTLD", "VCP"), 
       col=c("red","orange","green","blue","magenta"), cex=1.2) #add legend


#PCA on genes
pcaC9 <- prcomp(C9pat) #run pca
pcaCH <- prcomp(CHpat) #run pca
pcasALS <- prcomp(sALSpat) #run pca
pcaFTLD <- prcomp(FTLDpat) #run pca
pcaVCP <- prcomp(VCPpat) #run pca

dev.off() #removes previous plot
##
plot(pcaVCP$x[,1:2], pch=18, cex=1.25 , col="red",main = "PCA of Patient Microarray Datasets")

##
plot(pcaC9$x[,1:2], pch=18, cex=1.25 , col="red",main = "PCA of Patient Microarray Datasets", 
     xlim=c(-100000, 100000), ylim=c(-20000, 20000))
points(pcaCH$x[,1:2], pch=18, cex=1.25 , col="orange")
points(pcasALS$x[,1:2], pch=18, cex=1.25 , col="green")
points(pcaFTLD$x[,1:2], pch=18, cex=1.25 , col="blue")
points(pcaVCP$x[,1:2], pch=18, cex=1.25 , col="magenta")


library(rgl)
library(pca3d)

### C9orf72 ###
pca3d(pcaC9$x[,1:3]) #plot 3D
text3d(pcaC9$x[,1:3], text=rownames(pcaC9$x), adj=1.3, color="black", cex = 0.5) #add labels

### sALS ###
pcaSALS <- prcomp(t(exp_SALS.LCM))
pca3d(pcaSALS$x[,1:3])
text3d(pcaSALS$x[,1:3], text=rownames(pcaSALS$x), adj=1.3, color="black", cex = 0.5)


### CHMP2B ###
pcaCHMP2B <- prcomp(t(exp_CHMP2B.LCM))
pca3d(pcaCHMP2B$x[,1:3])
text3d(pcaCHMP2B$x[,1:3], text=rownames(pcaCHMP2B$x), adj=1.3, color="black", cex = 0.5)

### FTLD ###
pcaFTLD <- prcomp(t(FTLD_combo))
pca3d(pcaFTLD$x[,1:3])
text3d(pcaFTLD$x[,1:3], text=rownames(pcaFTLD$x), adj=1.3, color="black", cex = 0.5)

### VCP ###
pcaVCP <- prcomp(t(VCP))
pca3d(pcaVCP$x[,1:3])#plot 3D
text3d(pcaVCP$x[,1:3], text=rownames(pcaVCP$x), adj=1.3, color="black",cex = 0.5) #add labels
