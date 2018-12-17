load (file="C:/Users/sann1458/Dropbox/Transcriptome studies/Win260914/LCMdata.201015.R")

library (WGCNA)
library (pathprint)
library (sva)
library (biomaRt)
options(stringsAsFactors = FALSE)


static <- read.csv("C:/Users/sann1458/Dropbox/Win Network/Pathprint/staticFImodules_geneNames.csv", header=F)
enrich <- matrix (nrow=length(static[,1]), ncol=2)

for (i in 1: length(static[,1]))
{
  genes <- static[i,]
  enrich [i,1] <- as.character(genes [1])
  genes <- genes [-c(1)]
  blank <- which (genes=="")
  if (length (blank > 0))
  {
    genes <- genes [-c(blank)]
  }
  genes <- as.character(genes)
  
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mart_back <- getBM(attributes=c("hgnc_symbol", "affy_hg_u133_plus_2", "ensembl_transcript_id"), filters="hgnc_symbol", values=genes,  mart=mart)
  genes.affy <- unique (mart_back[,2])
  
  test <- length (which (Tracking.GM %in% genes.affy))
  
  ####check against random permutatuons
  
  m=10000
  r <- c(1:m)
  
  for (j in 1:m)
  {
    random <- sample (colnames(DATA), size=length (genes.affy), replace=F)
    random <- length (which (Tracking.GM %in% random))
    r[j] <- random
  }
  
  test1 <- which (r > test)
  enrich [i,2] <- (length(test1)/m)
  
}

write.table (enrich, file="C:/Users/sann1458/Dropbox/Win Network/191015/static.Tracking.GM.signif.txt")


####### if check in affy ids then only enrich BSG pathway and only one affy id when use tracking set SIGNIF ONLY (i.e. not GM)
######## FOR GM SIGNIF ONLY: ACTB, BSG, CREB1, CTSB, GRB2, RAP1A, RARA, SREBF1
