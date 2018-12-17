library(biomaRt)

NUG_list <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/TDP_fixed_NuggetGenes.txt")

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "ensembl_exon_id"), filters="hgnc_symbol", values=NUG_list,  mart=mart)

exonID <- mart_back$ensembl_exon_id

mart_back_Exon <- getBM(attributes =c("ensembl_exon_id","chromosome_name","exon_chrom_start","exon_chrom_end"), 
                        filters="ensembl_exon_id", values=exonID,  mart=mart)

Exon_Gene <- merge(mart_back, mart_back_Exon, by="ensembl_exon_id")

write.csv(Exon_Gene, "Claire_Network_Exons_chromnum.csv", row.names = F, quote = F)
