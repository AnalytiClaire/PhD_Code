## checking range of fold changes for 283 differentially expressed genes

DEGs <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_upanddown.txt")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
sals <- read.csv("sals_unique.csv")
ftld <- read.csv("ftld_unique.csv")
vcp <- read.csv("vcp_unique.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
rav <- read.csv("RAV_results_keepfiltering.csv")


geneplot <- list()
var <- vector()
sd <- vector()

for (i in 1:length(DEGs)){
  x <- DEGs[i]
  C9_gene <- C9[C9$Gene.Symbol %in% x,]
  sALS_gene <- sals[sals$Gene.Symbol %in% x,]
  FTLD_gene <- ftld[ftld$Gene.Symbol %in% x,]
  VCP_gene <- vcp[vcp$Gene.Symbol %in% x,]
  PET_gene <- pet[pet$hgnc_symbol %in% x,]
  RAV_gene <- rav[rav$hgnc_symbol %in% x,]
  
  foldchange <- cbind(C9_gene$Fold.Change, sALS_gene$Fold.Change, FTLD_gene$Fold.Change,
                      VCP_gene$Fold.Change, PET_gene$FoldChange, RAV_gene$FoldChange)
  colnames(foldchange) <- c("C9","sALS","FTLD","VCP","PET","RAV")
  foldchange <- t(foldchange)
  
  geneplot[[i]] <- foldchange
  var[[i]] <- var(foldchange)
  sd[[i]] <- sd(foldchange)
}

#####PLOT

library(ggplot2)
install.packages(c("ggplot2","RColorBrewer","scales"))
library(ggplot2); library(scales); library(grid); library(RColorBrewer)

boxplot(geneplot, 
        names = DEGs,
        las =2)

names(var) = DEGs
plot(sd,
     xlab = "Differentially Expressed Genes",
     ylab = "Standard Deviation")

SD = data.frame(row.names = DEGs,
            SD = sd)


# SD$over1 <- "FALSE"
# SD[c(56,58,106,129,132,149,233),2] <- "TRUE"
# SD$names <- rownames(SD)
# 
# for (i in 1:283){
#   if (SD$over1[i] == "FALSE"){
#     SD$names[i] <- "NA"
#   }
# }



sunflowerplot(var,
              xlab = "Differentially Expressed Genes",
              ylab = "Variance", 
              cex.lab = 2,
              cex = sd,
              mgp = c(2.5,0.7,0))
abline(0.5,0, col = "indianred1", lty = 5, lwd = 2)
abline(1.5,0, col = "indianred3", lty = 5, lwd = 2)
abline(2.5,0, col = "indianred4", lty = 5, lwd = 2)
text(var, DEGs, pos = 1, cex = sd)





boxplot(geneplot, ylim(c(4,-4)), main = "Fold Change of DEGs", 
        cex.main = 2, ylab = "Fold Change", cex.lab =1, 
        xlab = NULL)
abline(h=1, col = "indianred1", lty = 5, lwd = 2)
abline(h=1.5, col = "indianred3", lty = 5, lwd = 2)
abline(h=2, col = "indianred4", lty = 5, lwd = 2)

abline(h=-1, col = "steelblue1", lty = 5, lwd = 2)
abline(h=-1.5, col = "steelblue3", lty = 5, lwd = 2)
abline(h=-2, col = "steelblue4", lty = 5, lwd = 2)









