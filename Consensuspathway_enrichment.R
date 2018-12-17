# 
# library(pathprint)
# Pathprint_genes <- pathprint.Hs.gs
# 
# library (biomaRt)
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
# mart_attribute <- listAttributes(mart)
# 
# 
# pathprint_list <- list()
# for (i in 1:length(Pathprint_genes)){
#   annotation <- getBM(attributes=c("entrezgene", "hgnc_symbol"),
#                       filters = "entrezgene", values = Pathprint_genes[[i]], mart = mart)
#   pathprint_list[[i]] <- annotation$hgnc_symbol
# }
# 
# names(pathprint_list) <- names(Pathprint_genes)

# #0.8 up
# Z <- list(pathprint_list$`PPAR signaling pathway (KEGG)`, 
#           pathprint_list$`Glutamatergic synapse (KEGG)`,
#           pathprint_list$`Pancreatic secretion (KEGG)`,
#           pathprint_list$`Opioid Signalling (Reactome)`,
#           pathprint_list$`{CALM1,30} (Static Module)`)
# #0.8down
# Z <- list(pathprint_list$`Pentose phosphate pathway (KEGG)`, 
#           pathprint_list$`Other types of O-glycan biosynthesis (KEGG)`,
#           pathprint_list$`RNA polymerase (KEGG)`,
#           pathprint_list$`Proteasome (KEGG)`,
#           pathprint_list$`Androgen Receptor Signaling Pathway (Wikipathways)`, 
#           pathprint_list$`Regulatory RNA pathways (Reactome)`,
#           pathprint_list$`Signaling by Wnt (Reactome)`, 
#           pathprint_list$`IL-1 up reg. targets (Netpath)`, 
#           pathprint_list$`{EPRS,15} (Static Module)`,
#           pathprint_list$`{TFAP2A,25} (Static Module)`)
# #0.75 up
# Z <- list(pathprint_list$`ABC transporters (KEGG)`, 
#           pathprint_list$`PPAR signaling pathway (KEGG)`,
#           pathprint_list$`Aldosterone-regulated sodium reabsorption (KEGG)`,
#           pathprint_list$`Physiological and Pathological Hypertrophy  of the Heart (Wikipathways)`,
#           pathprint_list$`ACE Inhibitor Pathway (Wikipathways)`, 
#           pathprint_list$`Myogenesis (Reactome)` ,
#           pathprint_list$`Signaling by FGFR (Reactome)`, 
#           pathprint_list$`{CALM1, 30} (Static Module)`, 
#           pathprint_list$`{PAK2, 10} (Static Module)`)
# 
# #0.75 down
# Z <- list(pathprint_list$`Pentose phosphate pathway (KEGG)`, 
#           pathprint_list$`Other types of O-glycan biosynthesis (KEGG)`,
#           pathprint_list$`RNA polymerase (KEGG)`,
#           pathprint_list$`Proteasome (KEGG)`,
#           pathprint_list$`Androgen Receptor Signaling Pathway (Wikipathways)`, 
#           pathprint_list$`Regulatory RNA pathways (Reactome)`,
#           pathprint_list$`Signaling by Wnt (Reactome)`, 
#           pathprint_list$`{EPRS,15} (Static Module)`, 
#           pathprint_list$`{TFAP2A,25} (Static Module)`)
#



#down_path
# downpath <- list(pathprint_list$`{ABL1,15} (Static Module)`, 
#           pathprint_list$`{ACY1,11} (Static Module)`,
#           pathprint_list$`{ARRB2,743} (Static Module)`,
#           pathprint_list$`{CHRNA1,13} (Static Module)`,
#           pathprint_list$`{DVL1L1,17} (Static Module)`, 
#           pathprint_list$`{HRAS,27} (Static Module)`,
#           pathprint_list$`{TCF3,20} (Static Module)`, 
#           pathprint_list$`Aflatoxin B1 metabolism (Wikipathways)`, 
#           pathprint_list$`Alzheimer's disease (KEGG)`,
#           pathprint_list$`Blood Clotting Cascade (Wikipathways)`,
#           pathprint_list$`Collecting duct acid secretion (KEGG)`,
#           pathprint_list$`Diabetes pathways (Reactome)`,
#           pathprint_list$`Estrogen signaling pathway (Wikipathways)`,
#           pathprint_list$`Glucocorticoid &amp; Mineralcorticoid Metabolism (Wikipathways)`,
#           pathprint_list$`Influenza A (KEGG)`,
#           pathprint_list$`Phototransduction (KEGG)`,
#           pathprint_list$`Serotonin Receptor 4/6/7 and NR3C Signaling (Wikipathways)`,
#           pathprint_list$`Signaling by GPCR (Reactome)`,
#           pathprint_list$`Signaling by Insulin receptor (Reactome)`,
#           pathprint_list$`SNARE interactions in vesicular transport (KEGG)`,
#           pathprint_list$`Statin Pathway (Wikipathways)`,
#           pathprint_list$`Striated Muscle Contraction (Wikipathways)`,
#           pathprint_list$`Tyrosine metabolism (KEGG)`,
#           pathprint_list$`Urea cycle and metabolism of amino groups (Wikipathways)`)



#up_path
# uppath <- list(pathprint_list$`{BRCA1,28} (Static Module)`, 
#           pathprint_list$`{F2,46} (Static Module)`,
#           pathprint_list$`{FCGR2B,50} (Static Module)`,
#           pathprint_list$`{RB1,11} (Static Module)`,
#           pathprint_list$`{VCP,17} (Static Module)`, 
#           pathprint_list$`ABC transporters (KEGG)`,
#           pathprint_list$`Alanine, aspartate and glutamate metabolism (KEGG)`, 
#           pathprint_list$`Aminoacyl-tRNA biosynthesis (KEGG)`, 
#           pathprint_list$`Basal transcription factors (KEGG)`,
#           pathprint_list$`Codeine and morphine metabolism (Wikipathways)`,
#           pathprint_list$`Endometrial cancer (KEGG)`,
#           pathprint_list$`Fat digestion and absorption (KEGG)`,
#           pathprint_list$`Fatty Acid Biosynthesis (Wikipathways)`,
#           pathprint_list$`Glutathione metabolism (KEGG)`,
#           pathprint_list$`Glycine, serine and threonine metabolism (KEGG)`,
#           pathprint_list$`Id Signaling Pathway (Wikipathways)`,
#           pathprint_list$`IL-2 down reg. targets (Netpath)`,
#           pathprint_list$`Metabolic pathways (KEGG)`,
#           pathprint_list$`Mitochondrial LC-Fatty Acid Beta-Oxidation (Wikipathways)`,
#           pathprint_list$`Muscle contraction (Reactome)`,
#           pathprint_list$`Phagosome (KEGG)`,
#           pathprint_list$`Vitamin B12 Metabolism (Wikipathways)`,
#           pathprint_list$`Vitamin D synthesis (Wikipathways)`)

load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint_biomart.RData")

downpath <- subset(pathprint_list, names(pathprint_list) %in% down_pathways)
uppath <- subset(pathprint_list, names(pathprint_list) %in% up_pathways)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Consensus_fingerprint")
upnames <- read.table("upregulated_pathways.txt")
upnames <- upnames$V1

downnames <-read.table("downregulated_pathways.txt")
downnames <- downnames$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
up <- read.table("intersect_up_1.txt")
up <- up$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
down <- read.table("intersect_down_1.txt")
down <- down$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
upanddown <- read.table("upanddown.txt")
upanddown <- upanddown$V1

#Load file with all genes
library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]
sym.genes <- t(sym.genes)
allgenes <- sym.genes[!duplicated(sym.genes),]


# run script
pathwayEnrichment <- hyperPathway(
  genelist = W$Cirulli,
  geneset = uppath,
  Nchip = length(allgenes))
  
pathwayEnrichment$ID <- downnames
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Consensus_fingerprint")
write.csv(pathwayEnrichment, "downpath_upgene.csv")
  