
library(DESeq2)
library(biomaRt)
#### Analysis of Rav TDP-43 data ####

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data")
Counts <- read.csv("ravitsannotated_combined.counts.csv", row.names = 1)

#Assign condition
condition <- factor(c(rep("Condition1", 8), rep("Condition2", 13)))
condition

#Remove rows with all zeros
# Counts[rowSums(Counts) == 0,] <- NA
# Counts <- na.omit(Counts)
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
## Merge with raw and normalized count data
resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)


names(resdata_raw)[1] <- "Gene"
head(resdata_raw)
genes <- as.vector(resdata_raw$Gene)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)
result <- merge(resdata_raw, mart_back, by.x = "Gene", by.y = "ensembl_gene_id")
genesort <- result[order(result$pvalue),]
genesort<-subset(genesort, hgnc_symbol!="")


setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2")
write.csv(genesort, "RAV_results_GSM.csv")

# Sig.padj <- subset(genesort, subset=(padj < 0.05))
# Sig.padj <- Sig.padj$hgnc_symbol
# 
# write.table(Sig.padj, "sig.padj.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


######

#### Analysis of Pet TDP-43 data ####

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data")
Data <- read.csv("Pet.annotated_combined.counts.csv", row.names = 1)
Data <- Data[,1:27]


#Assign condition
# condition <- factor(c(rep("Condition1", 9), rep("Condition2", 7)))
condition <- factor(c(rep("Condition1", 9), rep("Condition2", 18)))
condition

#make column of gene IDs into rownames


#Remove rows with 
Counts <- Data
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
## Merge with normalized count data
resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)


names(resdata_raw)[1] <- "Gene"
head(resdata_raw)
genes <- as.vector(resdata_raw$Gene)
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)
result <- merge(resdata_raw, mart_back, by.x = "Gene", by.y = "ensembl_gene_id")
genesort <- result[order(result$pvalue),]
genesort<-subset(genesort, hgnc_symbol!="")


setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2")
write.csv(genesort, "PET_sALS_results_GSM.csv")

Sig.padj <- subset(genesort, subset=(padj < 0.05))
Sig.padj <- Sig.padj$hgnc_symbol

write.table(Sig.padj, "sig.padj.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
