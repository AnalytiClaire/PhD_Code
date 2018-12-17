##RNA-Seq Gene Expression Analysis using Limma##

library(limma)
library(edgeR)
library(biomaRt)

analysis.name<-"Pet" #Label analysis
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/Petrucelli/")
# Counts <- read.table(file = 'GSE67196_Petrucelli2015_ALS_genes.rawcount.txt', header = TRUE)
# 
# write.csv(x = Counts, file = "counts_petrucelli.csv")

Counts <- read.csv(file = "PET_Ens_Counts.csv", header = TRUE)
# Counts <- as.matrix(Counts[-28])
m <- as.matrix(Counts[,-1])
rownames(m) <- Counts[,1]
Counts <- m

Counts[rowSums(Counts) == 0,] <- NA
Counts <- na.omit(Counts)
Countnum <- Counts

#DGElist
dge <- DGEList(counts=Countnum)
dge <- calcNormFactors(dge)
colnames(dge)

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


genesort <- result[order(result$P.Value),]
genesort <- subset(genesort, subset=(hgnc_symbol !=""))

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")
write.csv(genesort, file=paste(analysis.name, "rankedresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)


####################################################
##RNA-Seq Gene Expression Analysis using Limma##

library(limma)
library(edgeR)
library(biomaRt)

analysis.name<-"Rav_naomit" #Label analysis
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/Ravits")

Counts <- read.csv(file = "ravitsannotated_combined.counts.csv", header = TRUE)
# Counts <- as.matrix(Counts[-28])
m <- as.matrix(Counts[,-1])
rownames(m) <- Counts[,1]
Counts <- m

Counts[rowSums(Counts) == 0,] <- NA
Counts <- na.omit(Counts)
Countnum <- Counts

#DGElist
dge <- DGEList(counts=Countnum)
dge <- calcNormFactors(dge)
colnames(dge)

#Design
Treat<-factor(rep(c("Control", "Patient"),c(8,13)), levels=c("Control", "Patient"))
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


#Count tables from bcbio have ensembl gene IDs. This must be annotated with HGNC symbols

#Download the HGNC symbols and gene IDs using a vector containing the IDs from results

genes <- as.vector(result[,1])
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)

#Merge the tables using ensembl ID
result <- merge(result, mart_back, by.x = "Row.names", by.y = "ensembl_gene_id")


#order by p value
genesort <- result[order(result$P.Value),]
#remove rows with no gene symbol
genesort <- subset(genesort, subset=(hgnc_symbol !=""))
#If duplicated, take highest value
genesort <- genesort[!duplicated(genesort$hgnc_symbol),]


setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")
write.csv(genesort, file=paste(analysis.name, "rankedresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)



