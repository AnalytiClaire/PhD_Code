library(ggplot2)
library(gridExtra)
library(rgl)
library(pca3d)
# 
# x <- read.csv("PAD_uniqueresult.csv")
# rownames(x) <- x$hgnc_symbol
# x <- x[9:26]
# write.csv(x, "PAD_expressiononly.csv")

### Quality Control ###
setwd("/users/claireblue/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
BER_OA <- read.csv("BER_OA_expressiononly.csv", row.names = 1)
BER_RA <- read.csv("BER_RA_expressiononly.csv", row.names = 1)
BRO_RA <- read.csv("BRO_RA_expressiononly.csv", row.names = 1)
JEN_OA <- read.csv("JEN_OA_expressiononly.csv", row.names = 1)
JEN_RA <- read.csv("JEN_RA_expressiononly.csv", row.names = 1)
VRI_OA <- read.csv("VRI_OA_expressiononly.csv", row.names = 1)
PAD_OA <- read.csv("PAD_expressiononly.csv", row.names = 1)

HOR <- read.csv("MID1ExpressionOnly.csv", row.names = 1)
MOU <- read.csv("MID2ExpressionOnly.csv", row.names = 1)
NAK <- read.csv("MOR.SNExpressionOnly.csv", row.names = 1)



##### BER_OA #####
# PCA 2D #
analysis.name <- "BER_OA"
data <- BER_OA
names <- colnames(BER_OA)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2.5, col=c(rep("blue",10), rep("red",10)),main = "BER_OA 2D PCA Plot: PC1 vs PC2") #plot disease
legend(25000, 30000,pch=18, legend=c("Disease", "Control"), col=c("red","blue"), cex=1.5) #add legend
# text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("blue",10), rep("red",10))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("blue",10),
                c(rep("red", 10))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.5,
        par(mar = c(10, 5, 1, 1)))

##### BER_RA #####
# PCA 2D #
analysis.name <- "BER_RA"
data <- BER_RA
names <- colnames(BER_RA)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2.5 , col=c(rep("blue",10), rep("red",10)),main = "BER_RA 2D PCA Plot: PC1 vs PC2") #plot disease
legend(40000, 38000,pch=18, legend=c("Disease", "Control"), col=c("red","blue"), cex=1.5) #add legend
# text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("blue",10), rep("red",10))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("blue",10),
                c(rep("red", 10))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.5,
        par(mar = c(10, 5, 1, 1)))


##### BRO_RA #####
# PCA 2D #
analysis.name <- "BRO_RA"
data <- BRO_RA
names <- colnames(BRO_RA)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2.5 , col=c(rep("blue",7), rep("red",16)),main = "BRO_RA 2D PCA Plot: PC1 vs PC2") #plot disease
legend(62000, 96000,pch=18, legend=c("Disease", "Control"), col=c("red","blue"), cex=1.5) #add legend
# text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("blue",7), rep("red",16))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("blue",7),
                c(rep("red", 16))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.7,
        par(mar = c(10, 5, 1, 1)))


##### JEN_OA #####
# PCA 2D #
analysis.name <- "JEN_OA"
data <- JEN_OA
names <- colnames(JEN_OA)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2.5 , col=c(rep("blue",10), rep("red",10)),main = "JEN_OA 2D PCA Plot: PC1 vs PC2") #plot disease
legend(15000, 53000,pch=18, legend=c("Disease", "Control"), col=c("red","blue"), cex=1.5) #add legend
# text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("blue",10), rep("red",10))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("blue",10),
                c(rep("red", 10))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.7,
        par(mar = c(10, 5, 1, 1)))


##### JEN_RA #####
# PCA 2D #
analysis.name <- "JEN_RA"
data <- JEN_RA
names <- colnames(JEN_RA)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2.5 , col=c(rep("blue",10), rep("red",13)),main = "JEN_RA 2D PCA Plot: PC1 vs PC2") #plot disease
legend(32000, 43000,pch=18, legend=c("Disease", "Control"), col=c("red","blue"), cex=1.5) #add legend
# text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("blue",10), rep("red",13))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("blue",10),
                c(rep("red", 13))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.7,
        par(mar = c(10, 5, 1, 1)))


##### VRI_OA #####
# PCA 2D #
analysis.name <- "VRI_OA"
data <- VRI_OA
names <- colnames(VRI_OA)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2.5 , col=c(rep("blue",7), rep("red",10)),main = "VRI_OA 2D PCA Plot: PC1 vs PC2") #plot disease
legend(80000, 75000,pch=18, legend=c("Disease", "Control"), col=c("red","blue"), cex=1.5) #add legend
# text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("blue",7), rep("red",10))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("blue",7),
                c(rep("red", 10))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.7,
        par(mar = c(10, 5, 1, 1)))


##### PAD_OA #####
# PCA 2D #
analysis.name <- "PAD_OA"
data <- PAD_OA
names <- colnames(PAD_OA)
dev.off() #removes previous plot
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2.5 , col=c(rep("blue",8), rep("red",10)),main = "PAD_OA 2D PCA Plot: PC1 vs PC2") #plot disease
legend(150000, 90000,pch=18, legend=c("Disease", "Control"), col=c("red","blue"), cex=1.5) #add legend
# text(pca_2D$x[,PC], labels = names,adj = 1.1, cex = 0.5) #label data points

#3d PCA
pca <- prcomp(t(data)) #conduct PCA
pca3d(pca$x[,1:3], col=c(rep("blue",8), rep("red",10))) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("blue",8),
                c(rep("red", 10))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.7,
        par(mar = c(10, 5, 1, 1)))
