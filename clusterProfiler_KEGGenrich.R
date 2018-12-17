library(DOSE)
library(biomaRt)
genes <- x
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("hgnc_symbol", "entrezgene"), filters="hgnc_symbol", values=genes,  mart=mart)



de <- mart_back$entrezgene

library(clusterProfiler)
kk <- enrichKEGG(de, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", 
                 qvalueCutoff=0.1)
head(summary(kk))
