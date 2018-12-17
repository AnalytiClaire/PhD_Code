### Heatmap for Benchmark Results
library(gplots)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/")
pvalue <- read.csv("adjpvalue_Heatmap.csv", )
geneshm <- read.csv("numberofgenes_Heatmap.csv")

rownames(pvalue) <- pvalue$GeneLists
pvalue[,1] <- NULL
pvalue <- as.matrix(pvalue)


rownames(geneshm) <- geneshm$GeneLists
geneshm[,1] <- NULL
geneshm <- as.matrix(geneshm)


### Heatmap for P Values
my_palette <- colorRampPalette(c("red","orange", "yellow", "green"))(n = 399)n

col_breaks = c(seq(0,0.00000005,length=100), # for red
               seq(0.000000051,0.009,length=100),  # for yellow
               seq(0.0091,0.05,length=100),
               seq(0.051,1,length=100))

heatmap.2(pvalue,
          cellnote = pvalue,    # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # only draw a row dendrogram
          Colv="NA",            #turn off column clustering
          Rowv = "NA")          # turn off row clustering
dev.off()               # close the PNG device


### Heatmap for Intersect Numbers


my_palette <- colorRampPalette(c("red","orange", "yellow", "green"))(n = 399)n

col_breaks = c(seq(0,0.00000005,length=100), # for red
               seq(0.000000051,0.009,length=100),  # for yellow
               seq(0.0091,0.05,length=100),
               seq(0.051,1,length=100))

heatmap.2(pvalue,
          cellnote = pvalue,    # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # only draw a row dendrogram
          Colv="NA",            #turn off column clustering
          Rowv = "NA")          # turn off row clustering
dev.off()               # close the PNG device


