library(ggplot2)
library(gridExtra)
library(rgl)
library(pca3d)

### Quality Control ###
setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
LEW <- read.csv("LEWExpressionOnly.csv", row.names = 1)
MID3 <- read.csv("MID3ExpressionOnly.csv", row.names = 1)
MID4 <- read.csv("MID4ExpressionOnly.csv", row.names = 1)
MOR.FC <- read.csv("MOR.FCExpressionOnly.csv", row.names = 1)
DIJ <- read.csv("DIJExpressionOnly.csv", row.names = 1)
FFR <- read.csv("FFRExpressionOnly.csv", row.names = 1)
MID1 <- read.csv("MID1ExpressionOnly.csv", row.names = 1)
MID2 <- read.csv("MID2ExpressionOnly.csv", row.names = 1)
MOR.SN <- read.csv("MOR.SNExpressionOnly.csv", row.names = 1)
BOT1 <- read.csv("BOT1ExpressionOnly.csv", row.names = 1)
BOT2 <- read.csv("BOT2ExpressionOnly.csv", row.names = 1)


##### LEW #####
# PCA 2D #
analysis.name <- "LEW"
data <- LEW
names <- colnames(LEW)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",4), rep("red",6)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(16000, 14000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",4), rep("red",6))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",4),
                c(rep("red", 6))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.7,
        par(mar = c(10, 5, 1, 1)))



##### MID3 #####
# PCA 2D #
analysis.name <- "MID3"
data <- MID3
names <- colnames(MID3)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",15), rep("red",14)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(16000, -20000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.3) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",15), rep("red",14))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",15),
                c(rep("red", 14))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(12, 5, 1, 1)))

##### MID4 #####

# PCA 2D #
analysis.name <- "MID4"
data <- MID4
names <- colnames(MID4)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(2,3)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",20), rep("red",15)),main = "2D PCA Plot: PC2 vs PC3") #plot disease
legend(20000, 25000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.3) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",20), rep("red",15))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",20),
                c(rep("red", 15))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(10, 5, 1, 1)))

##### MOR.FC #####
# PCA 2D #
analysis.name <- "MOR.FC"
data <- MOR.FC
names <- colnames(MOR.FC)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",3), rep("red",5)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(20000, 25000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.4) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",3), rep("red",5))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##

dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",3),
                c(rep("red", 5))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(10, 5, 1, 1)))

##### DIJ #####
# PCA 2D #
analysis.name <- "DIJ"
data <- DIJ
names <- colnames(DIJ)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",8), rep("red",15)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(20000, 18000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.4) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",8), rep("red",15))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",8),
                c(rep("red", 15))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(12, 5, 1, 1)))

##### FFR #####
# PCA 2D #
analysis.name <- "FFR"
data <- FFR
names <- colnames(FFR)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",9), rep("red",16)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(15000, 10000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.4) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",9), rep("red",16))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",9),
                c(rep("red", 16))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(10, 5, 1, 1)))

##### MID1 #####
# PCA 2D #
analysis.name <- "MID1"
data <- MID1
names <- colnames(MID1)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",8), rep("red",10)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(2000, 8000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.4) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",8), rep("red",10))) #plot 3D
# text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",8),
                c(rep("red", 10))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(14, 5, 1, 1)))

##### MID2 #####
# PCA 2D #
analysis.name <- "MID2"
data <- MID2
names <- colnames(MID2)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",18), rep("red",11)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(10000, 32000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.4) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",18), rep("red",11))) #plot 3D
# text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",18),
                c(rep("red", 11))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(10, 5, 1, 1)))


##### MOR.SN #####
# PCA 2D #
analysis.name <- "MOR.SN"
data <- MOR.SN
names <- colnames(MOR.SN)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",15), rep("red",24)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(-20000, 8000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.4) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",15), rep("red",24))) #plot 3D
# text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",15),
                c(rep("red", 24))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3, 
        par(mar = c(10, 5, 1, 1)))


##### BOT1 #####
# PCA 2D #
analysis.name <- "BOT1"
data <- BOT1
names <- colnames(BOT1)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",4), rep("red",2)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(-20000, 8000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.8) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",4), rep("red",2))) #plot 3D
# text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",4),
                c(rep("red", 2))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1.5,
        par(mar = c(10, 5, 1, 1)))

##### BOT2 #####
# PCA 2D #
analysis.name <- "BOT2"
data <- BOT2
names <- colnames(BOT2)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=1.25 , col=c(rep("green",5), rep("red",3)),main = "2D PCA Plot: PC1 vs PC2") #plot disease
legend(-20000, 8000,pch=18, legend=c("Disease", "Control"), col=c("red","green"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names,adj = 1, cex = 0.8) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("green",5), rep("red",3))) #plot 3D
# text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("green",5),
                c(rep("red", 3))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1.5,
        par(mar = c(10, 5, 1, 1)))
