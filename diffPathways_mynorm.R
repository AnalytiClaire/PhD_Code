##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####

library(pathprint)
data(list = c("chipframe", "genesets","pathprint.Hs.gs","platform.thresholds", "pluripotents.frame"))
library(pathprintGEOData)

options(stringsAsFactors = FALSE)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/NormalisedExpressionMatrices/") #set working directory to location of data

####C9_LCM ######
# setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/C9orf72_LCM") #set working directory to location of data

exp_C9.LCM <- read.csv ("C9eset.csv", header=TRUE, row.names = 1) #assign the .csv file to a variable, column headers are true
# row.names (exp_C9.LCM) <- exp_C9.LCM[,1] #specify that first column contains gene names
# exp_C9.LCM<- exp_C9.LCM[,2:12] #specify that all other columns are gene expression data

C9.LCM_pathprint <- exprs2fingerprint(exp_C9.LCM, platform = "GPL570", species="human", progressBar=T) #takes the gene expression values and converts into ternary score (-1,0,1) #platform = microarray platform GEO ID
vec.c9 <- c(0,0,0,1,1,1,1,1,1,1,1)


####sals_lcm###

# setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/FUS_SALS_LCM_CELfiles")
exp_SALS.LCM <- read.csv ("sALSeset.csv", header=TRUE, row.names = 1)
# row.names (exp_SALS.LCM) <- exp_SALS.LCM[,1]
# exp_SALS.LCM<- exp_SALS.LCM[,2:11]

SALS.LCM_pathprint <- exprs2fingerprint (exp_SALS.LCM, platform = "GPL570", species="human", progressBar=T)
vec.sals <- c(0,0,0,1,1,1,1,1,1,1)


####FTLD###
# 
# setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/FTD-U.brain")
exp_FTLD <- read.csv ("FTLDeset.csv", header=TRUE, row.names = 1)
# row.names (FTLD) <- FTLD[,1]
# FTLD <- FTLD[,2:25]

#GPL571 = Affymetrix Human Genome U113A 2.0 array
FTLD_pathprint <- exprs2fingerprint (exp_FTLD, platform = "GPL571", species="human", progressBar=T)
vec.FTLD <- c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)


####VCP###

# setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/VCP.myopathy")
exp_VCP <- read.csv ("VCPeset.csv", header=TRUE, row.names = 1)
# row.names (VCP) <- VCP[,1]
# VCP <- VCP[,2:11]

VCP_pathprint <- exprs2fingerprint (exp_VCP, platform = "GPL570", species="human", progressBar=T)
vec.vcp <- c(0,0,0,1,1,1,1,1,1,1)


##DiffPathways##

thres <- 0.1


threshpathC9 <- diffPathways(C9.LCM_pathprint, vec.c9, thres)
threshpathsals <- diffPathways(SALS.LCM_pathprint, vec.sals, thres)
threshpathftld <- diffPathways(FTLD_pathprint, vec.FTLD, thres)
threshpathvcp <- diffPathways(VCP_pathprint, vec.vcp, thres)



###INTERSECT###

overlap <- Reduce(intersect, list(threshpathC9, threshpathsals, threshpathftld, threshpathvcp)) #selects pathways that are present in all data sets listed
print(overlap)

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/")
write.table(overlap, "TDP_Diffpathways_mynorm.txt", quote = T, col.names = F, row.names = F)
# write.csv(VCP_pathprint, file = "VCPpathprint.csv")
# 
# C9 <- as.numeric(C9.LCM_pathprint['TGF beta receptor down reg. targets (Netpath)',])
# C9
# CH <- CHMP2B.lcm



#### consensus


