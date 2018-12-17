setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15")

values <- c("Gene.Symbol", "Fold.Change") #take values required

#Load all normalised results from DGEA
exprsC9 <- read.csv("C9rankeduniqueresult.csv")
exprsC9 <- exprsC9[values]
colnames(exprsC9)[2] <- "C9Fold.Change"

exprsCH <- read.csv("CHrankeduniqueresult.csv")
exprsCH <- exprsCH[values]
colnames(exprsCH)[2] <- "CHFold.Change"

exprssALS <- read.csv("sALSrankeduniqueresult.csv")
exprssALS <- exprssALS[values]
colnames(exprssALS)[2] <- "sALSFold.Change"

exprsFTLD <- read.csv("FTLDrankeduniqueresult.csv")
exprsFTLD <- exprsFTLD[values]
colnames(exprsFTLD)[2] <- "FTLDFold.Change"

exprsVCP <- read.csv("VCPrankeduniqueresult.csv")
exprsVCP <- exprsVCP[values]
colnames(exprsVCP)[2] <- "VCPFold.Change"

setwd(dir = "/Users/clairegreen/Desktop/")

#load the list containing genes of interest
Genelist <- read.csv("MADEGs.csv")

#Merge data from all data sets
Genelist <- merge(Genelist, exprsC9, by.x = "Gene", by.y = "Gene.Symbol")
Genelist <- merge(Genelist, exprsCH, by.x = "Gene", by.y = "Gene.Symbol")
Genelist <- merge(Genelist, exprssALS, by.x = "Gene", by.y = "Gene.Symbol")
Genelist <- merge(Genelist, exprsFTLD, by.x = "Gene", by.y = "Gene.Symbol")
Genelist <- merge(Genelist, exprsVCP, by.x = "Gene", by.y = "Gene.Symbol")

rownames(Genelist) <- Genelist[,1]
Genelist[,1] <- NULL

Genelistup <- Genelist
Genelistdown <- Genelist

#Show genes that are all upregulated
Genelistup[Genelistup < 0] <- NA
Genelistup <- na.omit(Genelistup)

#Show genes that are all downregulated
Genelistdown[Genelistdown > 0] <- NA
Genelistdown <- na.omit(Genelistdown)

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15")
write.csv(Genelist, file = "microarrayregulationdirection.csv")


#################################




