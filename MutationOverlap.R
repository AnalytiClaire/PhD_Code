##P2TDP Gene Overlap###
library(hgu133plus2.db)
library(pathprint)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/GeneMANIA/C9orf72")
C920 <- scan(file = "C9+20.txt", what = "list")
C950 <- scan(file = "C9+50.txt", what = "list")
C9100 <- scan(file = "C9+100.txt", what = "list")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/GeneMANIA/CHMP2B")
CH20 <- scan(file = "CHMP2B+20.txt", what = "list")
CH50 <- scan(file = "CHMP2B+50.txt", what = "list")
CH100 <- scan(file = "CHMP2B+100.txt", what = "list")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/GeneMANIA/GRN")
GRN20 <- scan(file = "GRN+20.txt", what = "list")
GRN50 <- scan(file = "GRN+50.txt", what = "list")
GRN100 <- scan(file = "GRN+100.txt", what = "list")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Pathways_to_TDP-43/GeneMANIA/VCP")
VCP20 <- scan(file = "VCP+20.txt", what = "list")
VCP50 <- scan(file = "VCP+50.txt", what = "list")
VCP100 <- scan(file = "VCP+100.txt", what = "list")

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
Y <- read.csv(file = "AllgenesNO.csv", na.strings = c("", "NA)"))
Y <- as.list(Y)
Y<- lapply(Y, function(x) x[!is.na(x)])

#Load file with all genes

sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)])
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

sym.genes <- t(sym.genes)

allgenes <- sym.genes[!duplicated(sym.genes),]

# run script
pathwayEnrichment <- hyperPathway(
  genelist = VCP100,
  geneset = Y,
  Nchip = length(allgenes)
)


overlap <- Reduce(intersect, list(VCP100, Y$X8000))
print(overlap)
