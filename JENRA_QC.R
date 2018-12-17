library(ggplot2)
library(gridExtra)
library(rgl)
library(pca3d)

### Quality Control ###
setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
BER.RA <- read.csv("BER_RA_filteredresult.csv", row.names = 1)
rownames(BER.RA) <- BER.RA$Gene.Symbol
BER.RA <- BER.RA[,19:38]

JEN.RA <- read.csv("JEN_RA_filteredresult.csv", row.names = 1)
rownames(JEN.RA) <- JEN.RA$Gene.Symbol
JEN.RA <- JEN.RA[,19:41]

BRO <- read.csv("BRO_RA_filteredresult.csv", row.names = 1)
rownames(BRO) <- BRO$Gene.Symbol
BRO <- BRO[,20:42]

WAL.RA <- read.csv("WAL_RA_UniqueGene_DESeq2.csv", row.names = 1)
rownames(WAL.RA) <- WAL.RA$Row.names
WAL.RA <- WAL.RA[,8:182]

DOL <- read.csv("DOL_filteredresult.csv", row.names = 1)
rownames(DOL) <- DOL$Gene.Symbol
DOL <- DOL[,19:28]

setwd("/users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PARA/GeneExpression/Fixed_JENRA/")
##### BER.RA #####
# PCA 2D #
analysis.name <- "BER.RA"
data <- BER.RA
names <- colnames(BER.RA)
p <- 10
c <- 10

pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("BER.RA_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-15000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
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
col = col=c(rep("salmon",10), rep("darkturquoise",10))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")



##### JEN.RA #####
# PCA 2D #
analysis.name <- "JEN.RA"
data <- JEN.RA
names <- colnames(JEN.RA)
p <- 13
c <- 10
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("JEN.RA_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-8000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
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
col = col=c(rep("salmon",10), rep("darkturquoise",5), 
            rep("salmon",2), rep("darkturquoise",5), "salmon")

names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")

##### BRO #####

# PCA 2D #
analysis.name <- "BRO"
data <- BRO
names <- colnames(BRO)
p <- 16
c <- 7
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("BRO_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-16000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
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
col = col=c(rep("darkturquoise",2), rep("salmon",5), 
            rep("darkturquoise",2), rep("salmon",7), "darkturquoise",
            "salmon", rep("darkturquoise",2), rep("salmon",3))

names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")

##### WAL.RA #####
# PCA 2D #
analysis.name <- "WAL.RA"
data <- WAL.RA[,c(1,3:125,127:167,169:175)]
names <- colnames(data)
p <- 150
c <- 27
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("WAL.RA_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-2000000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
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
x <- labels(clust)

write.table(x, "WAl.ralabels.txt", row.names = F, col.names = F, quote = F)
#Set right colours
col = readLines("WAl.ralabels.txt")
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 0.2)
plot(clust, xlab = "", sub = "")

##### DOL #####
# PCA 2D #
analysis.name <- "DOL"
data <- DOL
names <- colnames(DOL)
p <- 5
c <- 5
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("DOL_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-18000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
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
            "darkturquoise", "darkturquoise","darkturquoise","salmon","salmon")
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")

