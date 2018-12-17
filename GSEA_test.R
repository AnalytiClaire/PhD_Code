#### GSEA using GSEABase ####

library(GSEABase)
library(hgu133plus2cdf)
library(GO.db)
library(Biobase)

##C9orf72##

#set working directory to location of data
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/C9orf72_LCM") 

#assign the .csv file to a variable AS A MATRIX (not data.frame), column headers are true

exp_C9.LCM <- read.csv ("eset_NineP_150612_exprs.csv", header=TRUE)
C9_Mat <- as.matrix(read.csv("eset_NineP_150612_exprs.csv", header=TRUE))

#specify that first column contains gene names
row.names(exp_C9.LCM) <- exp_C9.LCM[,1]

#specify that all other columns are gene expression data
exp_C9.LCM<- exp_C9.LCM[,2:12]

x = phenoData(exp_C9.LCM)

#create a minimal ExpressionSet object using the ExpressionSet constructor
minimalSet <- ExpressionSet(assayData=exp_C9.LCM
                            phenoData)

egs <- GeneSet(minimalSet[1:54675,], setIdentifier = )
  
  
  