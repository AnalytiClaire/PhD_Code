##Convert information using biomart

#Important terms:
# Affymetrix Human Genome U133 Plus 2.0 array probe set = affy_hg_u133_plus_2
# Affymetrix Human Genome U133A 2 array probe set = affy_hg_u133a_2
# Ensembl ID = ensembl_gene_ensembl
# HGNC symbol = hgnc_symbol
# Entrez ID = entrezgene

overlap <- read.table("/Users/clairegreen/Desktop/My TDP-43 Genes.txt", header = FALSE)
overlap <- overlap$V1

library (biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

analysis.name <- "overlap"

genes <- overlap
# genes <- genes$Probe.Set.ID
mart_attribute <- listAttributes(mart)
mart_filter <- listFilters(mart)
annotation <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "hgnc_symbol", values = genes, mart = mart)


annotation
setwd("/Users/clairegreen/Desktop/")
write.csv(annotation, file = "overlap_ens2hgnc.csv")


# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/GoodData/")
# write.csv(annotation, paste(analysis.name,"geneIDs.csv"), sep = "")

SNPS <- genes$SNP

snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "www.ensembl.org")
SNP_attribute <- listAttributes(snp.db)
SNP_filter <- listFilters(snp.db)
nt.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start", "chrom_strand","associated_gene",  "ensembl_gene_stable_id"),  
                    filters="ensembl_gene",   values=SNPS,    mart=snp.db)
