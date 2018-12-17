library(ggplot2)
library(gridExtra)
library(rgl)
library(pca3d)
library(dendroextras)

### Quality Control ###
setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/")
C9 <- read.csv("C9uniquegene_samples.csv", row.names = 1)
sals <- read.csv("sALSuniquegene_samples.csv", row.names = 1)
ftld <- read.csv("FTLDuniquegene_samples.csv", row.names = 1)
vcp <- read.csv("VCPuniquegene_samples.csv", row.names = 1)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_sALS_results_GSM.csv")
rownames(pet) <- pet$hgnc_symbol
pet <- pet[,9:35]
rav <- read.csv("RAV_results_GSM.csv")
rownames(rav) <- rav$hgnc_symbol
rav <- rav[,9:29]


setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/GeneExpression/")

dev.off() #removes previous plot
##### C9 #####
analysis.name <- "C9"
data <- C9
names <- colnames(C9)
p <- 8
c <- 3

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("C9_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-1000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off() #removes previous plot


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
col = col=c("salmon","salmon","darkturquoise","salmon","darkturquoise",
            "darkturquoise", "salmon","salmon","salmon","salmon","salmon")
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### sALS #####
analysis.name <- "sALS"
data <- sals
names <- colnames(sals)
p <- 7
c <- 3

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("sALS_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-4000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.7) #label data points
dev.off() #removes previous plot


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
col = col=c("salmon","salmon","salmon","salmon","salmon",
            "salmon", "darkturquoise","darkturquoise","darkturquoise","salmon")
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(2,0,0,0), cex = 1)
plot(clust, xlab = "", sub = "")
                 

##### FTLD #####
analysis.name <- "FTLD"
data <- ftld
names <- colnames(ftld)
p <- 16
c <- 8

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("FTLD_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-6000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.4) #label data points
dev.off() #removes previous plot

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
par(oma=c(2,0,0,0), cex = 0.6)
plot(clust, xlab = "", sub = "")

#Set right colours
col = col=c(rep("salmon",8),rep("darkturquoise",2),
            rep("salmon",2), rep("darkturquoise",2),
            rep("salmon",1), rep("darkturquoise", 2),
            rep("salmon",5), rep("darkturquoise", 2))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(2,0,0,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### VCP #####
analysis.name <- "VCP"
data <- vcp
names <- colnames(vcp)
p <- 7
c <- 3

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("VCP_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-8000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off() #removes previous plot


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
col = col=c(rep("salmon",5),rep("darkturquoise",3),
            rep("salmon",2))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(2,0,0,0), cex = 1)
plot(clust, xlab = "", sub = "")


##### PET #####
analysis.name <- "PET"
data <- pet
names <- colnames(pet)
p <- 27
c <- 9


# PCA 2D #
dev.off()
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("PET_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-900000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.4) #label data points
dev.off() #removes previous plot


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
col = col=c(rep("salmon",1), rep("darkturquoise",2),
            rep("salmon",1), rep("darkturquoise",1),
            rep("salmon",1), rep("darkturquoise",4), 
            rep("salmon",15), rep("darkturquoise",2))
names(col) <- names

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,1,0), cex = 0.9)
plot(clust, xlab = "", sub = "")


##### RAV ##### SAME FOR RAV
analysis.name <- "RAV"
data <- rav
names <- colnames(rav)
p <- 13
c <- 8

# PCA 2D #
pca_2D <- prcomp(t(data))
PC <- c(1,2)
pdf("RAV_PCA.pdf", width=11, height=8) 
plot(pca_2D$x[,PC], pch=18, cex=3, col=c(rep("darkturquoise",c), rep("salmon",p)),main = paste(analysis.name)) #plot disease
legend(max(pca_2D$x[,1])-800000, max(pca_2D$x[,2])-10,pch=18, legend=c("Disease", "Control"), col=c("salmon","darkturquoise"), cex=1) #add legend
text(pca_2D$x[,PC], labels = names, cex = 0.6) #label data points
dev.off() #removes previous plot

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
col = col=c("salmon", "salmon", "salmon", "darkturquoise", "darkturquoise",
            "salmon", "darkturquoise", "darkturquoise", "salmon",
            "darkturquoise", "salmon", "salmon", "salmon", "salmon", 
            "salmon", "salmon", "darkturquoise",  "salmon","darkturquoise",
            "salmon")

#Replot
clust=colour_clusters(hclust(dist(t(data))), col = col, k = NULL, h = TRUE)
set_leaf_colours(d=clust, col=col, col_to_set = "label")
par(oma=c(0,0,2,0), cex = 1)
plot(clust, xlab = "", sub = "")
