
### 2D Principle Component Analysis ###

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/Ravits/")
pet <- read.csv("counts_petrucelli.csv")
rav <- read.csv("GSE76220_ALS_LCM_RPKM.csv")

library(ggplot2)
library(gridExtra)

### C9orf72 ###
dev.off() #removes previous plot
C9d <- exp_C9.LCM[,c(1:8)] #take expression data from disease columns
C9c <- exp_C9.LCM[,c(9:11)] #take expression data from control columns
pcaC9d <- prcomp(t(C9d)) #run pca
pcaC9c <- prcomp(t(C9c))

p_C9d <- plot(pcaC9d$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), 
              ylim=c(-200, 350),main = "C9orf72") #plot disease
p_C9c <- points(pcaC9c$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), 
              ylim=c(-200, 350)) #plot control
legend(260, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.7) #add legend

text(pcaC9d$x[,1:2], labels = rownames(pcaC9d$x),adj = 1.1, cex = 0.5) #label data points
text(pcaC9c$x[,1:2], labels = rownames(pcaC9c$x),adj = 1.1, cex = 0.5)



# ### CHMP2B ###
# dev.off()
# CHMP2Bd <- exp_CHMP2B.LCM[,c(1:7)]
# CHMP2Bc <- exp_CHMP2B.LCM[,c(8:10)]
# pca_CHMP2Bd <- prcomp(t(CHMP2Bd))
# pca_CHMP2Bc <- prcomp(t(CHMP2Bc))
# 
# p_CHMP2Bd <- plot(pca_CHMP2Bd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), 
#                   ylim=c(-200, 350), main = "CHMP2B")
# p_CHMP2Bc <- points(pca_CHMP2Bc$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), 
#                   ylim=c(-200, 350))
# legend(260, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.7)
# 
# text(pca_CHMP2Bd$x[,1:2], labels = rownames(pca_CHMP2Bd$x),adj = 1.1, cex = 0.5)
# text(pca_CHMP2Bc$x[,1:2], labels = rownames(pca_CHMP2Bc$x),adj = 1.1, cex = 0.5)
# 


### sALS ###
dev.off()
sALSd <- exp_SALS.LCM[,c(1:3)]
sALSc <- exp_SALS.LCM[,c(4:10)]
pcasALSd <- prcomp(t(sALSd))
pcasALSc <- prcomp(t(sALSc))

p_sALSd <- plot(pcasALSd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), 
                ylim=c(-200, 350), main = "sALS")
p_sALSc <- points(pcasALSc$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), 
                ylim=c(-200, 350))
legend(260, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.7)

text(pcasALSd$x[,1:2], labels = rownames(pcasALSd$x),adj = 1.1, cex = 0.5)
text(pcasALSc$x[,1:2], labels = rownames(pcasALSc$x),adj = 1.1, cex = 0.5)



### FTLD ###
dev.off()
pca_FTLD_c <- prcomp(t(FTLD.cont))
pca_FTLD_prgn <- prcomp(t(FTLD.prgn))
pca_FTLD_sftd <- prcomp(t(FTLD.sftd))

p_FTLD_prgn <- plot(pca_FTLD_prgn$x[,1:2], pch=18, cex=1.25 , col="orange", xlim=c(-300, 400), 
                    ylim=c(-200, 350),main = "FTLD")
p_FTLD_sftd <- points(pca_FTLD_sftd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), 
                    ylim=c(-200, 350))
p_FTLD_c <- points(pca_FTLD_c$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), 
                    ylim=c(-200, 350))
legend(260, 340,pch=18, legend=c("PRGN", "SFTD", "Control"), col=c("orange","red","green"), cex=0.7)

text(pca_FTLD_c$x[,1:2], labels = rownames(pca_FTLD_c$x),adj = 1.1, cex = 0.5)
text(pca_FTLD_prgn$x[,1:2], labels = rownames(pca_FTLD_prgn$x),adj = 1.1, cex = 0.5)
text(pca_FTLD_sftd$x[,1:2], labels = rownames(pca_FTLD_sftd$x),adj = 1.1, cex = 0.5)



### VCP ###
dev.off()
VCPd <- VCP[,c(1:3)]
VCPc <- VCP[,c(4:10)]
pca_VCPd <- prcomp(t(VCPd))
pca_VCPc <- prcomp(t(VCPc))

p_VCPd <- plot(pca_VCPd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), 
               ylim=c(-200, 350),main = "VCP")
p_VCPc <- points(pca_VCPc$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), 
               ylim=c(-200, 350))
legend(260, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.7)

text(pca_VCPd$x[,1:2], labels = rownames(pca_VCPd$x),adj = 1.1, cex = 0.5)
text(pca_VCPc$x[,1:2], labels = rownames(pca_VCPc$x),adj = 1.1, cex = 0.5)



### EL ###
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/EL_TDP-43/")
EL <- read.csv("EL_results_15052017_normalised2.csv")
dev.off() #removes previous plot
ELd <- EL[,c(9:15)] #take expression data from disease columns
ELc <- EL[,c(16:22)] #take expression data from control columns
pcaC9d <- prcomp(t(C9d)) #run pca
pcaC9c <- prcomp(t(C9c))

p_C9d <- plot(pcaC9d$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), 
              ylim=c(-200, 350),main = "C9orf72") #plot disease
p_C9c <- points(pcaC9c$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), 
                ylim=c(-200, 350)) #plot control
legend(260, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.7) #add legend

text(pcaC9d$x[,1:2], labels = rownames(pcaC9d$x),adj = 1.1, cex = 0.5) #label data points
text(pcaC9c$x[,1:2], labels = rownames(pcaC9c$x),adj = 1.1, cex = 0.5)

#### COMBINE IN DIFFERENT PANES ####
dev.off()
par(mfrow=c(3,2))
p1 <- plot(pcaC9d$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), ylim=c(-200, 350),
                    main = "C9orf72")
points(pcaC9c$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), ylim=c(-200, 350))
legend(250, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.5)


p2 <- plot(pca_CHMP2Bd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), ylim=c(-200, 350),
                        main = "CHMP2B")
points(pca_CHMP2Bc$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), ylim=c(-200, 350))
legend(250, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.5)

p3 <- plot(pcasALSd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), ylim=c(-200, 350),
           main = "sALS")
points(pcasALSc$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), ylim=c(-200, 350))
legend(250, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.5)

p4 <- plot(pca_FTLD_prgn$x[,1:2], pch=18, cex=1.25 , col="orange", xlim=c(-300, 400), ylim=c(-200, 350),
           main = "FTLD")
points(pca_FTLD_sftd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_FTLD_c$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), ylim=c(-200, 350))
legend(250, 340,pch=18, legend=c("PRGN", "SFTD", "Control"), col=c("orange","red","green"), cex=0.5)

p5 <- p_VCPd <- plot(pca_VCPd$x[,1:2], pch=18, cex=1.25 , col="red", xlim=c(-300, 400), ylim=c(-200, 350),
                     main = "VCP")
points(pca_VCPc$x[,1:2], pch=18, cex=1.25,col="green", xlim=c(-300, 400), ylim=c(-200, 350))
legend(250, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.5)


#### COMBINE IN ONE PANE ####

###DISEASE VS CONTROL####
dev.off()

plot(pcaC9d$x[,1:2], pch=18, cex=1.25, col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pcaC9c$x[,1:2], pch=18, cex=1.25, col="green", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_CHMP2Bd$x[,1:2], pch=18, cex=1.25, col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_CHMP2Bc$x[,1:2], pch=18, cex=1.25, col="green", xlim=c(-300, 400), ylim=c(-200, 350))
points(pcasALSd$x[,1:2], pch=18, cex=1.25, col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pcasALSc$x[,1:2], pch=18, cex=1.25, col="green", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_FTLD_prgn$x[,1:2], pch=18, cex=1.25, col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_FTLD_sftd$x[,1:2], pch=18, cex=1.25, col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_FTLD_c$x[,1:2], pch=18, cex=1.25, col="green", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_VCPd$x[,1:2], pch=18, cex=1.25, col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pca_VCPc$x[,1:2], pch=18, cex=1.25, col="green", xlim=c(-300, 400), ylim=c(-200, 350))

legend(250, 340,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=0.5)

#### ALL ####
dev.off()

pcaC9 <- prcomp(t(exp_C9.LCM))
pcaSALS <- prcomp(t(exp_SALS.LCM))
pcaCHMP2B <- prcomp(t(exp_CHMP2B.LCM))
pcaFTLD <- prcomp(t(FTLD_combo))
pcaVCP <- prcomp(t(VCP))

dev.off()
plot(pcaC9$x[,1:2], pch=18, col="red", xlim=c(-300, 400), ylim=c(-200, 350))
points(pcaSALS$x[,1:2],pch=18, col="blue")
points(pcaCHMP2B$x[,1:2],pch=18, col="green")
points(pcaFTLD$x[,1:2],pch=18, col="magenta")
points(pcaVCP$x[,1:2],pch=18, col="cyan")
legend(260, 340,pch=18, legend=c("C9orf72","sALS","CHMP2B","FTLD","VCP"), col=c("red", "blue", "green", "red", "cyan"), cex=0.7)

