##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####

library (pathprint)
options(stringsAsFactors = FALSE)


####C9_LCM ######
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/C9orf72_LCM") #set working directory to location of data
exp_C9.LCM <- read.csv ("eset_NineP_150612_exprs1.csv", header=TRUE) #assign the .csv file to a variable, column headers are true
row.names (exp_C9.LCM) <- exp_C9.LCM[,1] #specify that first column contains gene names
exp_C9.LCM<- exp_C9.LCM[,2:12] #specify that all other columns are gene expression data

C9.LCM_pathprint <- exprs2fingerprint(exp_C9.LCM, platform = "GPL570", species="human", progressBar=T) #takes the gene expression values and converts into ternary score (-1,0,1) #platform = microarray platform GEO ID
vec.c9 <- c(1,1,1,1,1,1,1,1,0,0,0)


####CHMP2B_LCM ######
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/CHMP2B")
exp_CHMP2B.LCM <- read.csv ("eset_CHMP2B_250615_exprs nooutlier.csv", header=TRUE)
row.names (exp_CHMP2B.LCM) <- exp_CHMP2B.LCM[,1]
exp_CHMP2B.LCM<- exp_CHMP2B.LCM[,2:10]

CHMP2B.LCM_pathprint <- exprs2fingerprint (exp_CHMP2B.LCM, platform = "GPL570", species="human", progressBar=T)
vec.ch <- c(1,1,1,0,0,0,0,0,0)


####sals_lcm###

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/FUS_SALS_LCM_CELfiles")
exp_SALS.LCM <- read.csv ("eset_SALS_LCM_260615_exprs.csv", header=TRUE)
row.names (exp_SALS.LCM) <- exp_SALS.LCM[,1]
exp_SALS.LCM<- exp_SALS.LCM[,2:11]

SALS.LCM_pathprint <- exprs2fingerprint (exp_SALS.LCM, platform = "GPL570", species="human", progressBar=T)
vec.sals <- c(0,0,0,1,1,1,1,1,1,1)


####FTLD###

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/FTD-U.brain")
FTLD <- read.csv ("FTLD_expr_tdp43.csv", header=TRUE)
row.names (FTLD) <- FTLD[,1]
FTLD <- FTLD[,2:25]

#GPL571 = Affymetrix Human Genome U113A 2.0 array
FTLD_pathprint <- exprs2fingerprint (FTLD, platform = "GPL571", species="human", progressBar=T)
vec.FTLD <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)


####VCP###

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/VCP.myopathy")
VCP <- read.csv ("eset_VCP.myopathy_170715_exprs.csv", header=TRUE)
row.names (VCP) <- VCP[,1]
VCP <- VCP[,2:11]

VCP_pathprint <- exprs2fingerprint (VCP, platform = "GPL570", species="human", progressBar=T)
vec.vcp <- c(0,0,0,1,1,1,1,1,1,1)


##DiffPathways##

thres <- 0.16


c9.lcm <- diffPathways(C9.LCM_pathprint, vec.c9, thres)
CHMP2B.lcm <- diffPathways(CHMP2B.LCM_pathprint, vec.ch, thres)
SALS.lcm <- diffPathways(SALS.LCM_pathprint, vec.sals, thres)
FTLD_FCx <- diffPathways(FTLD_pathprint, vec.FTLD, thres)
VCP.m <- diffPathways(VCP_pathprint, vec.vcp, thres)


### LRRK2/PARKIN Analysis ###




###INTERSECT###

overlap <- Reduce(intersect, list(c9.lcm, CHMP2B.lcm, SALS.lcm, FTLD_FCx, VCP.m)) #selects pathways that are present in all data sets listed
print(overlap)

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/")

write.csv(VCP_pathprint, file = "VCPpathprint.csv")

C9 <- as.numeric(C9.LCM_pathprint['TGF beta receptor down reg. targets (Netpath)',])
C9
CH <- CHMP2B.lcm



