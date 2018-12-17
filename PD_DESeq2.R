#DESeq2 for PD RNA-Seq
library(DESeq2)
library(biomaRt)
options(scipen=999)
# #Input from salmon
# library(tkWidgets)
# 
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/bcbio/salmon.quants/")
# 
# #run program to choose .CEL files from directory
# temp = list.files(pattern="*.csv")
# myfiles = lapply(temp, read.csv)
# 
# salmon.counts <- matrix(nrow = 194993, ncol = 73)
# rownames(salmon.counts) <- myfiles[[1]]$Name
# 
# for (i in 1:length(myfiles)){
#   salmon.counts[,i] <- myfiles[[i]]$NumReads
# }
# 
# # names <- unlist(strsplit(temp, "quant.csv"))
# # write.table(names, "GSMnumbers.txt", quote = F, row.names = F, col.names = F)
# 
# names <- readLines("GSMnumbers.txt")
# colnames(salmon.counts) <- names
# 
# write.csv(salmon.counts, "salmon.counts.csv")

#####################################


Counts <- read.csv("salmon.counts.int.csv", row.names = 1, stringsAsFactors = F)
trans <- rownames(Counts)


#Assign condition
condition <- factor(c(rep("Condition1", 44), rep("Condition2", 29)))
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
## Merge with raw and normalized count data
# resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)


names(resdata_norm)[1] <- "Gene"
head(resdata_norm)
genes <- as.vector(resdata_norm$Gene)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_transcript_id", "hgnc_symbol"), filters="ensembl_transcript_id", values=genes,  mart=mart)
result <- merge(resdata_norm, mart_back, by.x = "Gene", by.y = "ensembl_transcript_id")
genesort <- result[order(result$pvalue),]
genesort<-subset(genesort, hgnc_symbol!="")
uniqueresult <- genesort[!duplicated(genesort[,81]),]


setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(uniqueresult, "DUM_UniqueGene_DESeq2.csv")
