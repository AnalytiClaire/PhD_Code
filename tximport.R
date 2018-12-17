### Salmon Processing ##

library(tximport)
library(biomaRt)
setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Data/WAL")
files <- list.files()

ENST <- readLines("../ENSTtoHGNC.txt")
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("ensembl_transcript_id","hgnc_symbol"), filters="ensembl_transcript_id", values=ENST,  mart=mart)


txi.salmon <- tximport(files=files, type = "salmon", tx2gene = mart_back)
head(txi.salmon$counts)

counts <- txi.salmon$counts
colnames(counts) <- files

setwd("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Data/")
write.csv(counts, "WAL_counts.csv")
