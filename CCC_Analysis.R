##FINDING TOP GENES IN CCC PATHWAY FROM INDIVIDUAL MUTATIONS###

############################ Mutation DEGS #####################################
setwd(dir = "/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")

C9_array <- read.csv("C9rankeduniqueresult.csv")
C9_gene <- subset(C9_array, P.Value < 0.05)
C9_sym <- C9_gene$Gene.Symbol

#C9_1000 <- C9_gene[1:1000, drop=TRUE]

CH <- read.csv("CHrankeduniqueresult.csv")
CH_gene <- subset(CH, P.Value < 0.05)
CH_sym <- CH_gene$Gene.Symbol

#CH_gene <- CH$Gene.Symbol
#CH_1000 <- CH_gene[1:1000, drop=TRUE]

FTLD <- read.csv("GRNrankeduniqueresult.csv")
FTLD_gene <- subset(FTLD, P.Value < 0.05)
FTLD_sym <- FTLD_gene$Gene.Symbol

#FTLD_gene <- FTLD$Gene.Symbol
#FTLD_1000 <- FTLD_gene[1:1000, drop=TRUE]

VCP <- read.csv("VCPrankeduniqueresult.csv")
VCP_gene <- subset(VCP, P.Value < 0.05)
VCP_sym <- VCP_gene$Gene.Symbol

VCP_gene <- VCP$Gene.Symbol
VCP_2000 <- VCP_gene[1:4000, drop=TRUE]


############################ TDP-43 interaction genes #####################################
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
Q <- read.table(file = "Pasterkamp_TDP43.txt")
q <- Q$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
R <- read.table(file = "Taylor_TDP43.txt")
r <- R$V1

############################ PCxN Gene Lists #####################################
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results")

amoebiasis <- read.table(file = "Amoebiasis_genes.txt")
amo <- amoebiasis$V1

cytokine <- read.table(file = "Cytokine-cytokine_receptor_interaction_genes.txt")
cyto <- cytokine$V1 

IL4down <- read.table(file = "IL-4_down_reg._targets_genes.txt")
IL4 <- IL4down$V1
#VIM

integrin <- read.table(file = "Integrin_cell_surface_interactions_genes.txt")
int <- integrin$V1

malaria <- read.table(file = "Malaria_genes.txt")
mal <- malaria$V1

adipogenesis <- read.table(file = "Adipogenesis_genes.txt")
adi <- adipogenesis$V1

IL6up <- read.table(file = "IL-6_up_reg._targets.txt")
IL6 <- IL6up$V1

phagosome <- read.table(file = "Phagosome_genes.txt")
phag <- phagosome$V1
#"TUBA1B" "TUBB" 

selenium <- read.table(file = "Selenium_pathway_genes.txt")
sel <- selenium$V1
#ALB

vitaminb12 <- read.table(file = "Vitamin_B12_metabolism_genes.txt")
vit <- vitaminb12$V1
#ALB

############################ Complement and Coagulation Cascade genes #####################################
CCC <- read.table("/users/clairegreen/Desktop/Desktop Folder/CCC_KEGG_genes.txt")
CCC <- CCC$V1

############################ All pathway genes ###############################
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint 25.04.16/hyperPathway/All.Pathways/PathprintPathways(29)/")

Z <- read.csv(file = "pathprintgenes.csv", na.strings = c("", "NA)"))
Z <- as.list(Z)
Z <- lapply(Z, function(x) x[!is.na(x)])

############################ Common DEGs #####################################
DEGs <- read.table("/users/clairegreen/Desktop/Desktop Folder/TDP-43 DEGs.txt")
DEGs <- DEGs$V1


########################## DEGs from Guillaume's TDP43 cells ################
G_DEGs <- read.csv("/users/clairegreen/Desktop/Guillaume_TDP43_DEGs.csv", na.strings = c("", "NA)"))
G_DEGs<- as.list(G_DEGs)
G_DEGs <- lapply(G_DEGs, function(x) x[!is.na(x)])




############################ Overlap #####################################
overlap <- Reduce(intersect, list(e, CCC))
overlap


######## Hypergeometric Test #####

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
  genelist = C9_sym,
  geneset = Z,
  Nchip = length(allgenes)
)
setwd (dir = "/Users/clairegreen/Desktop/")
write.csv(pathwayEnrichment, file = "VCPpathways_enrich.csv")



x <- setdiff(diff1, diff2)
print(x)
