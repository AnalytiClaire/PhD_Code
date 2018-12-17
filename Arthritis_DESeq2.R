#DESeq2 for Arthritis RNA-Seq
library(DESeq2)
library(biomaRt)

### PAD ###
setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Data/PAD")
Counts <- read.csv("PAD_Counts.int.csv", row.names = 1, stringsAsFactors = F)
Counts[,12] <- NULL #Remove GSM3130542 as it did not pass quality control
trans <- rownames(Counts)


#Assign condition - controls/untreated first, patients/treated second
condition <- factor(c(rep("Condition1", 8), rep("Condition2", 9)))
condition

#Keep rows with at least 3 scores of 10 or above
keep <- rowSums(Counts>=10) >= 3
Counts<-Counts[keep,]

#Create a coldata frame
coldata <- data.frame(row.names=colnames(Counts), condition)
coldata

#Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData=Counts, colData=coldata, design=~condition)
dds

#Run DEseq2 pipleline
dds <- DESeq(dds)

res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$pvalue), ]
## Merge with raw or normalized count data
# resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

#bcbio comes with ensembl IDs so find the HGNC symbols using biomart
names(resdata_norm)[1] <- "Gene" #rename gene name column
head(resdata_norm)
genes <- as.vector(resdata_norm$Gene)
#set mart to use ensembl
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org") 
#Ask for ensembl and hgnc symbol matches using the id you have (ensembl ID) as a filter, and your list as the values
mart_back <- getBM(attributes =c("ensembl_transcript_id", "hgnc_symbol"), filters="ensembl_transcript_id", values=genes,  mart=mart)
#merge result with your result table
result <- merge(resdata_norm, mart_back, by.x = "Gene", by.y = "ensembl_transcript_id")
result$"Fold.Change"<-2^result$log2FoldChange 
result$"Fold.Change"[result$"Fold.Change"<1]<-(-1)/result$"Fold.Change"[result$"Fold.Change"<1] #converts log fold change into a linear value above or below 0
#order by pvalue
genesort <- result[order(result$pvalue),]
#remove rows with no hgnc symbol
genesort<-subset(genesort, hgnc_symbol!="")
#remove duplicate values (highest DE is kept)
uniqueresult <- genesort[!duplicated(genesort[,25]),]


setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression")
write.csv(uniqueresult, "PAD_uniqueresult.csv")



### Arthritis DESeq2 ###

#DESeq2 for PD RNA-Seq
library(DESeq2)
library(biomaRt)
options(scipen=999)
# #Input from salmon

setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Data/")
Counts <- read.csv("WAL_counts.csv", stringsAsFactors = F)
Counts <- Counts[!duplicated(Counts[,1]),]
rownames(Counts) <- Counts[,1]
Counts[,1] <- NULL
trans <- rownames(Counts)


#### WAL_OA ###
#Assign condition

WAL_OA <- Counts[,1:50]

condition <- factor(c(rep("Condition1", 28), rep("Condition2", 22)))
condition

#Keep rows with at least 3 scores of 10 or above
keep <- rowSums(WAL_OA>=10) >= 3
WAL_OA<-WAL_OA[keep,]

#Create a coldata frame
coldata <- data.frame(row.names=colnames(WAL_OA), condition)
coldata

#Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData=WAL_OA, colData=coldata, design=~condition)
dds

#Run DEseq2 pipleline
dds <- DESeq(dds)

res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$pvalue), ]
## Merge with raw and normalized count data
# resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)


# names(resdata_norm)[1] <- "Gene"
# head(resdata_norm)
# genes <- as.vector(resdata_norm$Gene)
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
# mart_back <- getBM(attributes =c("ensembl_transcript_id", "hgnc_symbol"), filters="ensembl_transcript_id", values=genes,  mart=mart)
result <- resdata_norm
result$"Fold.Change"<-2^result$log2FoldChange 
result$"Fold.Change"[result$"Fold.Change"<1]<-(-1)/result$"Fold.Change"[result$"Fold.Change"<1]
genesort <- result[order(result$pvalue),]
genesort<-subset(genesort, Row.names!="")
uniqueresult <- genesort[!duplicated(genesort[,1]),]


setwd(dir = "/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
write.csv(uniqueresult, "WAL_OA_UniqueGene_DESeq2.csv")


#### WAL_RA ###
#Assign condition

WAL_RA <- Counts[,c(1:28, 51:197)]

condition <- factor(c(rep("Condition1", 28), rep("Condition2", 147)))
condition

#Keep rows with at least 3 scores of 10 or above
keep <- rowSums(WAL_RA>=10) >= 3
WAL_RA<-WAL_RA[keep,]

#Create a coldata frame
coldata <- data.frame(row.names=colnames(WAL_RA), condition)
coldata

#Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData=WAL_RA, colData=coldata, design=~condition)
dds

#Run DEseq2 pipleline
dds <- DESeq(dds)

res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$pvalue), ]
## Merge with raw and normalized count data
# resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)


# names(resdata_norm)[1] <- "Gene"
# head(resdata_norm)
# genes <- as.vector(resdata_norm$Gene)
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
# mart_back <- getBM(attributes =c("ensembl_transcript_id", "hgnc_symbol"), filters="ensembl_transcript_id", values=genes,  mart=mart)
result <- resdata_norm
result$"Fold.Change"<-2^result$log2FoldChange 
result$"Fold.Change"[result$"Fold.Change"<1]<-(-1)/result$"Fold.Change"[result$"Fold.Change"<1]
genesort <- result[order(result$pvalue),]
genesort<-subset(genesort, Row.names!="")
uniqueresult <- genesort[!duplicated(genesort[,1]),]


setwd(dir = "/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/GeneExpression/")
write.csv(uniqueresult, "WAL_RA_UniqueGene_DESeq2.csv")

