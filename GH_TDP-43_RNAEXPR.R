##RNA-Seq Gene Expression Analysis using Limma##

library(limma)
library(edgeR)
library(biomaRt)
library(plyr)

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h")
# Counts <- read.table(file = 'GSE67196_Petrucelli2015_ALS_genes.rawcount.txt', header = TRUE)
# 
# write.csv(x = Counts, file = "counts_petrucelli.csv")

Counts <- read.csv(file = "combined.counts.csv", header = TRUE)

# Counts[Counts == 0] <- NA
# # Counts[Counts<30] <- NA
# Counts <- na.omit(Counts)
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


#####
analysis.name<-"TRL" #Label analysis

#Select experiment
Countnum <- TRL

Countnum[Countnum == 0] <- NA
Countnum <- na.omit(Countnum)
# rownames(Countnum)<-Countnum[,1]
# Countnum[,1] <- NULL

#DGElist
dge <- DGEList(counts=Countnum)
dge <- calcNormFactors(dge)

#Design
Treat<-factor(rep(c("Control", "Patient"),c(3,3)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(Countnum)
design

#Voom transformation
v <- voom(dge,design,plot=FALSE)

#Limma fitting
fit <- lmFit(v,design)
fit <- eBayes(fit)
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(Countnum)) #"BH" adjust for multiple hypothesis testing
result <- merge(result, Countnum, by="row.names", all=TRUE)
#result <- result[,1:7]

#Count tables from bcbio have ensembl gene IDs. This must be annotated with HGNC symbols

#Download the HGNC symbols and gene IDs using a vector containing the IDs from results

genes <- as.vector(result[,1])
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)

# library(org.Hs.eg.db)
# library(GeneNetworkBuilder)

#Merge the tables using ensembl ID
result <- merge(result, mart_back, by.x = "Row.names", by.y = "ensembl_gene_id")
# result[,1] <- NULL




#### Take median value for gene duplicates ###########
result3 <- ddply(result,"hgnc_symbol", numcolwise(median, (result$adj.P.Val)))
#result3 <- aggregate(result, by=list("Gene.Symbol"), FUN=median)

genesort <- result3[order(result3$adj.P.Val),]

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/results_2017_02_15_(moved filter)/")
write.csv(genesort, file=paste(analysis.name, "_GH_HEK_48hr.csv", sep=""), row.names=FALSE, quote = FALSE)









# 
# uniqueresult <- result[!duplicated(result$hgnc_symbol),]
# rownames(uniqueresult) <- uniqueresult$hgnc_symbol
# genesort <- uniqueresult[order(uniqueresult$adj.P.Val),]
# 
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
# write.csv(genesort, file=paste(analysis.name, "EXPRSrankeduniqueresult.csv", sep=""), sep="\t", row.names=TRUE, quote = FALSE)

# topgene <- genesort[1:1000,]
# write.csv(x = topgene, file = paste(analysis.name,"_ap_1000.csv", sep = ""))
# topgene <- genesort[1:2000,]
# write.csv(x = topgene, file = paste(analysis.name,"_ap_2000.csv", sep = ""))
# topgene <- genesort[1:3000,]
# write.csv(x = topgene, file = paste(analysis.name,"_ap_3000.csv", sep = ""))
# topgene <- genesort[1:4000,]
# write.csv(x = topgene, file = paste(analysis.name,"_ap_4000.csv", sep = ""))
# topgene <- genesort[1:5000,]
# write.csv(x = topgene, file = paste(analysis.name,"_ap_5000.csv", sep = ""))
