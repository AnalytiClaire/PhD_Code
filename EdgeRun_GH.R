### Edge Run

library(edgeR)
library(data.table)
library(edgeRun)
library(biomaRt)
library(plyr)

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h")

Counts <- read.csv(file = "combined.counts.csv", header = TRUE)

rownames(Counts)<-Counts[,1]
Counts[,1] <- NULL

## Divide into different analyses
TRL <- Counts[,1:6]
WCT <- Counts[,13:18]
CYT <- Counts[,19:24]
GFP_L <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","GRASPS.TRL.GFPLow.1",
                   "GRASPS.TRL.GFPLow.2","GRASPS.TRL.GFPLow.3")]
GFP_H <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","GRASPS.TRL.GFPHigh.1",
                   "GRASPS.TRL.GFPHigh.2","GRASPS.TRL.GFPHigh.3")]
GRASPS_VS_WCT <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","WCT.Control.1",
                           "WCT.Control.2","WCT.Control.3")]
GRASPS_VS_CYT <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","CytTrc.Control.1",
                           "CytTrc.Control.2","CytTrc.Control.3")]
WCT_VS_CYT <- Counts[,c("WCT.Control.1","WCT.Control.2","WCT.Control.3","CytTrc.Control.1",
                        "CytTrc.Control.2","CytTrc.Control.3")]
#Select experiment
analysis.name <- "CYT"
countdata <- CYT

samples <- c(1,1,1,2,2,2)

list <- DGEList(counts = countdata, lib.size = colSums(countdata), samples = samples, remove.zeros = TRUE)
d <- calcNormFactors(list)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

de.edgeR <- exactTest(d)
table <- as.data.frame(topTags(de.edgeR, n = length(de.edgeR$table[,1])))
table[,5] <- 0
colnames(table)[5] <- "Ens_ID"
table$Ens_ID <- rownames(table)

genes <- as.vector(rownames(table))
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)
result <- merge(table, mart_back, by.x = "Ens_ID", by.y = "ensembl_gene_id")
result_med <- ddply(result,"hgnc_symbol", numcolwise(median, (result$PValue)))
#result3 <- aggregate(result, by=list("Gene.Symbol"), FUN=median)

genesort <- result[order(result$PValue),]

## Write results
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/EdgeR")
write.csv(genesort, file=paste(analysis.name, "_diffexpr-results.csv", sep=""), row.names=FALSE, quote = FALSE)

CG_DEGs <- genesort

### 

CG_sig_EdgeR <- subset(CG_DEGs, subset=(FDR < 0.05))
CG_TRL_gene_EdgeR <- CG_sig_EdgeR$hgnc_symbol
write.table(CG_TRL_gene_EdgeR, file=paste(analysis.name, "_diffexpr-genes.txt", sep=""), row.names=FALSE, quote = FALSE)

overlap <- Reduce(intersect, list())
print(overlap)

