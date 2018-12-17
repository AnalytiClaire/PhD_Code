#### Analysis of EL's TDP-43 data ####

# TDPneg = nuclei absent of TDP-43
# TDPpos = nuclei containing TDP-43
options(scipen=999)

library(DESeq2)
library(biomaRt)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/EL_TDP-43/")
load("Count_Matrix_Genes.rda")

EnsNames <- rownames(counts.g)
EnsNames <- EnsNames[grepl("^ENS", EnsNames)]

counts.g.ens <- subset(counts.g, rownames(counts.g) %in% EnsNames)
counts.g.hgnc <- subset(counts.g, !(rownames(counts.g) %in% EnsNames))


ens_genes <- row.names(counts.g.ens)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=ens_genes,  mart=mart)


hgnc_genes <- row.names(counts.g.hgnc)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back_hgnc <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="hgnc_symbol", values=hgnc_genes,  mart=mart)


counts_ens <- merge(counts.g, mart_back, by.x = 0, by.y = "ensembl_gene_id")
counts_hgnc <- merge(counts.g, mart_back_hgnc, by.x = 0, by.y = "hgnc_symbol")


counts_ens <- counts_ens[,c(16,1,3,5,7,9,11,13,15,2,4,6,8,10,12,14)]
counts_hgnc <- counts_hgnc[,c(1,16,3,5,7,9,11,13,15,2,4,6,8,10,12,14)]

colnames(counts_ens)[1] <- "hgnc_symbol"
colnames(counts_hgnc)[1] <- "hgnc_symbol"
colnames(counts_ens)[2] <- "Ensembl_ID"
colnames(counts_hgnc)[2] <- "Ensembl_ID"

counts_all <- rbind(counts_ens, counts_hgnc)
counts_all<- subset(counts_all, subset=(hgnc_symbol !="")) #if no gene symbol, discount


exp_info <- data.frame(condition = factor(c(rep("1", 7), rep("2", 7))),
                       patientID = factor(c(1,2,3,4,5,6,7,1,2,3,4,5,6,7)),
                       Sex = factor(c("F","M","F","F","M","M","M")), 
                       Disease = factor(c("FTLD","FTLD-ALS","FTLD","FTLD","FTLD","FTLD-ALS","FTLD-ALS")))

  §
#Rows must have at least 3 samples with scores of 10 or higher 
keep <- rowSums(counts_all>=10) >= 3
counts_all<-counts_all[keep,]

counts_all_data <- counts_all
rownames(counts_all_data) <- counts_all$Ensembl_ID
counts_all_data <- counts_all_data[,3:16]

#Create a coldata frame
coldata <- data.frame(row.names=colnames(counts_all_data), exp_info)
coldata

#Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData=counts_all_data, colData=coldata, design=~patientID + conditions)
dds

#Run DEseq2 pipleline
dds <- DESeq(dds)

res <- results(dds)
table(res$padj<0.05)
## Order by p-value
res <- res[order(res$pvalue), ]

# ## Merge with raw count data
# resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)

## Merge with normalised count data
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
# resdata_norm_data <- resdata_norm[,8:21]
genenames <- counts_all[,1:2]
result <- merge(resdata_norm, genenames, by.x = "Row.names", by.y = "Ensembl_ID")
result <- result[,c(1,22,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
colnames(result)[1] <- "Ensembl_ID"

genesort <- result[order(result$pvalue),]
genesort <- genesort[!duplicated(genesort$hgnc_symbol),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/EL_TDP-43/2017_05_31/")
write.csv(genesort, "EL_results_31052017", row.names = F)

Sig.padj <- subset(genesort, subset=(padj < 0.05))
Sig.padj.gene <- Sig.padj$hgnc_symbol
Sig.padj.gene <- Sig.padj.gene[!duplicated(Sig.padj.gene)]

Sig.padj.up <- subset(Sig.padj, subset=(log2FoldChange > 0))
Sig.padj.up.gene <- Sig.padj.up$hgnc_symbol
Sig.padj.up.gene <- Sig.padj.up.gene[!duplicated(Sig.padj.up.gene)]

Sig.padj.down <- subset(Sig.padj, subset=(log2FoldChange < 0))
Sig.padj.down.gene <- Sig.padj.down$hgnc_symbol
Sig.padj.down.gene <- Sig.padj.down.gene[!duplicated(Sig.padj.down.gene)]

Sig.padj.both <- c(Sig.padj.down.gene, Sig.padj.up.gene)

write.table(Sig.padj.both, "sig_padj_genenames.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Sig.padj.up.gene, "sig_padj_up_genenames.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Sig.padj.down.gene, "sig_padj_down_genenames.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)




