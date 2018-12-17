## 3D Principle Component Analysis ##

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/Quality Control Analysis/")
load("TDP-43 Analysis Environment.RData")


#DISCLAIMER
#To use rgl you have to have X11 installed. If you get the error that X11 is not found (usually
#mac), then go to XQuartx and download it. After that, load the rgl source("y") again. 

library(rgl)
library(pca3d)

### C9orf72 ###
pcaC9 <- prcomp(t(exp_C9.LCM)) #conduct PCA
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
