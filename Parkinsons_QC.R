library(ggplot2)
library(gridExtra)
library(rgl)
library(pca3d)
library(dendroextras)

### Quality Control ###
setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
LEW <- read.csv("LEWExpressionOnly.csv", row.names = 1)
FFR <- read.csv("FFRExpressionOnly.csv", row.names = 1)
DIJ <- read.csv("DIJExpressionOnly.csv", row.names = 1)
MID1 <- read.csv("MID1ExpressionOnly.csv", row.names = 1)
MID2 <- read.csv("MID2ExpressionOnly.csv", row.names = 1)
MID3 <- read.csv("MID3ExpressionOnly.csv", row.names = 1)
MID4 <- read.csv("MID4ExpressionOnly.csv", row.names = 1)
MOR.FC <- read.csv("MOR.FCExpressionOnly.csv", row.names = 1)
MOR.SN <- read.csv("MOR.SNExpressionOnly.csv", row.names = 1)
DUM <- read.csv("DUM_UniqueGene_DESeq2.csv")
rownames(DUM) <- DUM$hgnc_symbol
DUM <- DUM[,9:81]
BOT1 <- read.csv("BOT1ExpressionOnly.csv", row.names = 1)
BOT2 <- read.csv("BOT2ExpressionOnly.csv", row.names = 1)

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/QC/")

dev.off() #removes previous plot
##### LEW #####
analysis.name <- "LEW"
data <- LEW
names <- colnames(LEW)
p <- 6
c <- 4

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("LEW_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-3000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.8) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
library(dendroextras)
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(2,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c("salmon","salmon","salmon","darkturquoise","darkturquoise",
            "salmon", "darkturquoise","salmon","darkturquoise","salmon")
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")
########################
##### C_GSM488113 ######
########################


##### FFR #####
analysis.name <- "FFR"
data <- FFR
names <- colnames(FFR)
p <- 16
c <- 9

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("FFR_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-4000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(2,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c("salmon", rep("darkturquoise", 4), "salmon", rep("darkturquoise", 2),
            rep("salmon",9), rep("darkturquoise",3), rep("salmon", 5))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(2,0,0,0), cex = 1)
plot(clust, xlab = "", sub = "")

########################
##### P_GSM184366 ######
########################


##### DIJ #####
analysis.name <- "DIJ"
data <- DIJ
names <- colnames(DIJ)
p <- 15
c <- 8

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("DIJ_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-8000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(2,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("salmon",3),rep("darkturquoise",7),
            rep("salmon",5), rep("darkturquoise",1),
            rep("salmon",7))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### MID1 #####
analysis.name <- "MID1"
data <- MID1
names <- colnames(MID1)
p <- 10
c <- 8

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("MID1_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-1500, max(pca_2D$x[,2]-100),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1.2,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("salmon",2),rep("darkturquoise",1),
            rep("salmon",2), rep("darkturquoise",1),
            rep("salmon",5), rep("darkturquoise",1), 
            rep("salmon",1),rep("darkturquoise",5))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,0,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### MID2 #####
analysis.name <- "MID2"
data <- MID2
names <- colnames(MID2)
p <- 11
c <- 18


# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("MID2_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-14000, max(pca_2D$x[,2]+800),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression",
        las = 2,
        cex = 0.3,
        cex.axis = 0.8,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")


#Set right colours
col = col=c(rep("salmon",1), rep("darkturquoise",6),
            rep("salmon",1), rep("darkturquoise",1),
            rep("salmon",4), rep("darkturquoise",8), 
            rep("salmon",4), rep("darkturquoise",1),
            rep("salmon",1), rep("darkturquoise",2))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.9)
plot(clust, xlab = "", sub = "")

########################
##### P_GSM508732 ######
########################

##### MID3 #####
analysis.name <- "MID3"
data <- MID3
names <- colnames(MID3)
p <- 14
c <- 15

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("MID3_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-6000, max(pca_2D$x[,2]+1000),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("salmon", 5), rep("darkturquoise", 2),
            rep("salmon", 1), rep("darkturquoise", 1),
            rep("salmon", 1), rep("darkturquoise", 1),
            rep("salmon", 1), rep("darkturquoise", 5),
            rep("salmon", 4), rep("darkturquoise", 4),
            rep("salmon", 2), rep("darkturquoise", 2))

names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")



##### MID4 #####
analysis.name <- "MID4"
data <- MID4
names <- colnames(MID4)
p <- 15
c <- 20

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("MID4_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-6000, max(pca_2D$x[,2]+1000),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 2), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 2),
            rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 3), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 3),
            rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 3), rep("salmon", 1),
            rep("darkturquoise", 3))

names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### MOR.FC #####
analysis.name <- "MOR.FC"
data <- MOR.FC
names <- colnames(MOR.FC)
p <- 5
c <- 3

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("MOR.FC_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-4000, max(pca_2D$x[,2]),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("salmon", 3), rep("darkturquoise", 2),
            rep("salmon", 1), rep("darkturquoise", 1),
            rep("salmon", 1))

names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")





##### MOR.SN #####
analysis.name <- "MOR.SN"
data <- MOR.SN
names <- colnames(MOR.SN)
p <- 24
c <- 15

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("MOR.SN_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-30000, max(pca_2D$x[,2]),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 0.8,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("darkturquoise", 2), rep("salmon", 10),
            rep("darkturquoise", 3), rep("salmon", 2),
            rep("darkturquoise", 2), rep("salmon", 12),
            rep("darkturquoise", 8))


names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 0.8)
plot(clust, xlab = "", sub = "")

##### DUM #####
analysis.name <- "DUM"
data <- DUM
names <- colnames(DUM)
p <- 29
c <- 44

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("DUM_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-500000, max(pca_2D$x[,2]),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.2,
        cex.axis = 0.8,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = readLines("alsihf.txt")


names(col) <- labels(clust)

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

### 711 ###

##### BOT1 #####
analysis.name <- "BOT1"
data <- BOT1
names <- colnames(BOT1)
p <- 2
c <- 4

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("BOT1_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(6000,1850,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("darkturquoise", 3), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 1))


names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### BOT2#####
analysis.name <- "BOT2"
data <- BOT2
names <- colnames(BOT2)
p <- 3
c <- 5

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("BOT2_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-2000, max(pca_2D$x[,2]),pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off()


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(0,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("darkturquoise", 1), rep("salmon", 1),
            rep("darkturquoise", 3), rep("salmon", 1),
            rep("darkturquoise", 1), rep("salmon", 1))


names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")


####################################################
setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

AMA <- read.csv("AMAfilteredresult.csv")
rownames(AMA) <- AMA$Gene.Symbol
AMA <- AMA[,21:372]

RON <- read.csv("RONfilteredresult.csv")
rownames(RON) <- RON$Gene.Symbol
RON <- RON[,20:78]

ATP13A2 <- read.csv("AMA_ATP13A2filteredresult.csv")
rownames(ATP13A2) <- ATP13A2$Gene.Symbol
ATP13A2 <- ATP13A2[,21:208]

PRKN <- read.csv("AMA_PARKINfilteredresult.csv")
rownames(PRKN) <- PRKN$Gene.Symbol
PRKN <- PRKN[,21:216]

PINK1 <- read.csv("AMA_PINK1filteredresult.csv")
rownames(PINK1) <- PINK1$Gene.Symbol
PINK1 <- PINK1[,21:215]

dev.off() #removes previous plot
##### AMA #####
analysis.name <- "AMA"
data <- AMA
names <- colnames(AMA)
p <- 169
c <- 183

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-8000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
library(dendroextras)
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(2,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = readLines("AMACol.txt")
            
            
            
            
names(col) <- labels(clust)

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, sub = "")



##### RON #####
analysis.name <- "RON"
data <- RON
names <- colnames(RON)
p <- 40
c <- 19

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-6000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.8) #label data points


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1.2,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(2,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = readLines("RonCol.txt")
names(col) <- labels(clust)

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(2,0,0,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### AMA_ATP13A2 #####
analysis.name <- "ATP13A2"
data <- ATP13A2
names <- colnames(ATP13A2)
p <- 5
c <- 183

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=3 , col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-6000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.4) #label data points


## Box plots ##
dev.off()
logdata <- log2(data)
boxplot(logdata,
        col = c(rep("darkturquoise",c),
                c(rep("salmon", p))),
        names = names,
        ylab = "Log2 Expression", 
        las = 2, 
        cex = 0.3,
        cex.axis = 1.2,
        par(mar = c(10, 5, 1, 1)))

# Cluster #
#Find out order of samples
clust=colour_clusters(hclust(dist(t(data))), col = c(rep("black", ncol(data))), k = NULL, h = TRUE)
set_leaf_colours(d=clust, col='black', col_to_set = "label")
par(oma=c(2,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = readLines("ATP13A2.txt")
names(col) <- labels(clust)

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(2,0,0,0), cex = 1)
plot(clust, xlab = "", sub = "")
