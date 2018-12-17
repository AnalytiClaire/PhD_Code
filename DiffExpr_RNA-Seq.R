##RNA-Seq Gene Expression Analysis using Limma##

library(limma)
library(edgeR)
library(biomaRt)

analysis.name<-"GH_HEK_48h" #Label analysis
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/Petrucelli/")
# Counts <- read.table(file = 'GSE67196_Petrucelli2015_ALS_genes.rawcount.txt', header = TRUE)
# 
# write.csv(x = Counts, file = "counts_petrucelli.csv")

Counts <- read.csv(file = "PETasMatrix.csv", header = TRUE)
m <- as.matrix(Counts[,-1])
rownames(m) <- Counts[,1]
Counts <- m

Counts <- as.matrix(Counts)
rownames(Counts)<-Counts[,1]
Counts <- Counts[,2:27]
# Counts[,1] <- NULL

Counts[rowSums(Counts) == 0,] <- NA
# Counts[Counts<30] <- NA
Counts <- na.omit(Counts)



# Counts<-subset(Counts, subset=(GeneID !="NA")) #if no gene symbol, discount

# Countszero <-subset(Counts, subset=(row !=0))
# Countszero <- apply(Counts, 1, function(row) all(row !="NA"))
# Counts <- Counts[Countszero,]



#remove gene symbols from end column
# Countnum <- Counts[,1:26]
# Counts <- data.matrix(Counts)

# Countnum <- read.csv(file = "pet.counts.clean.csv")

#DGElist
dge <- DGEList(counts=Countnum)
dge <- calcNormFactors(dge)

#Design
Treat<-factor(rep(c("Control", "Patient"),c(9,17)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(Countnum)
design

#Voom transformation
v <- voom(dge,design,plot=FALSE)

#Limma fitting
fit <- lmFit(v,design)
fit <- eBayes(fit)
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(Countnum)) #"BH" adjust for multiple hypothesis testing
result <- merge(result, Counts, by="row.names", all=TRUE)
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

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/MedianGenes")
write.csv(genesort, file=paste(analysis.name, "HGNCrankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)












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
