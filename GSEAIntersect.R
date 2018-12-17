##GSEA intersect##

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/C9orf72 without outlier/")
C9_GSEA <- read.csv(file = "gsea_report_for_PAT_1448468512104.csv")
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/CHMP2B without oulier/")
CHMP2B_GSEA <- read.csv(file = "gsea_report_for_PAT_1448622432572.csv")
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/sALS/")
sALS_GSEA <- read.csv(file = "gsea_report_for_PAT_1448533319830.csv")
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/FTLD/")
FTLD_GSEA <- read.csv(file = "gsea_report_for_PAT_1448633043887.csv")
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/VCP/")
VCP_GSEA <- read.csv(file = "gsea_report_for_PAT_1448542310559.csv")\


C9l <- C9_GSEA[,1]
CHMP2Bl <- CHMP2B_GSEA[,1]
sALSl <- sALS_GSEA[,1]
FTLDl <- FTLD_GSEA[,1]
VCPl <- VCP_GSEA[,1]

overlap <- Reduce(intersect, list(C9l, CHMP2Bl, sALSl, FTLDl, VCPl))
print(overlap)

