#### Sample Clustering to Detect Outliars using WGCNA ####

library(WGCNA)


### C9orf72 ###
# Display the current working directory
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/C9orf72_LCM")
options(stringsAsFactors = FALSE);

#Read in data set
exp_C9.LCM <- read.csv ("eset_NineP_150612_exprs.csv", header=TRUE)

row.names (exp_C9.LCM) <- exp_C9.LCM[,1] #specify that first column contains gene names
exp_C9.LCM<- exp_C9.LCM[,2:12] #specify that all other columns are gene expression data
exp_C9.LCM <- t(exp_C9.LCM) #transpose data set so that samples are rows and genes are columns (required by sampleTree)

###Check that there are any excessive missing values and identification of outlier microarray###
gsg = goodSamplesGenes(exp_C9.LCM, verbose = 3);
gsg$allOK

#If this comes out as False, look up documentation here http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf

###Cluster the samples###

C9_sampleTree = hclust(dist(exp_C9.LCM), method = "average");


### CHMP2B ###
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/CHMP2B/")
options(stringsAsFactors = FALSE);

exp_CHMP2B <- read.csv ("eset_CHMP2B_250615_exprs.csv", header=TRUE)

row.names (exp_CHMP2B) <- exp_CHMP2B[,1] 
exp_CHMP2B<- exp_CHMP2B[,2:11] 
exp_CHMP2B <- t(exp_CHMP2B)

gsg = goodSamplesGenes(exp_CHMP2B, verbose = 3);
gsg$allOK


###Cluster the samples###

CHMP2B_sampleTree = hclust(dist(exp_CHMP2B), method = "average");


### sALS ###
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/FUS_SALS_LCM_CELfiles/")
options(stringsAsFactors = FALSE);

exp_sALS <- read.csv ("eset_SALS_LCM_260615_exprs.csv", header=TRUE)

row.names (exp_sALS) <- exp_sALS[,1] 
exp_sALS<- exp_sALS[,2:11] 
exp_sALS <- t(exp_sALS) 

gsg = goodSamplesGenes(exp_sALS, verbose = 3);
gsg$allOK

###Cluster the samples###

sALS_sampleTree = hclust(dist(exp_sALS), method = "average");


### FTLD ###
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/FTD-U.brain/")
options(stringsAsFactors = FALSE);

FTLD <- read.csv ("eset_FTD.U.brain_170715_exprs.csv", header=TRUE)

row.names (FTLD) <- FTLD[,1]
FTLD <- FTLD[,2:57]
FTLD.cont <- FTLD[,c(1,4,5,7,8,11,13,15)]
FTLD.prgn <- FTLD[,c(18,20,22,24,27,30)]
FTLD.sftd <- FTLD[,c(33,35,38,41,44,45,48,50,52,55)]
FTLD_combo <- cbind(FTLD.cont, FTLD.prgn, FTLD.sftd) #combine control, progranulin and sporadic samples

FTLD_combo <- t(FTLD_combo) 

gsg = goodSamplesGenes(FTLD_combo, verbose = 3);
gsg$allOK

###Cluster the samples###

FTLD_sampleTree = hclust(dist(FTLD_combo), method = "average");


### VCP ###
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/VCP.myopathy/")
options(stringsAsFactors = FALSE);

VCP <- read.csv ("eset_VCP.Myopathy_170715_exprs.csv", header=TRUE)

row.names (VCP) <- VCP[,1] 
VCP<- VCP[,2:11] 
VCP <- t(VCP) 

gsg = goodSamplesGenes(VCP, verbose = 3);
gsg$allOK

###Cluster the samples###

VCP_sampleTree = hclust(dist(VCP), method = "average");


### PLOT CLUSTER DIAGRAMS ###
dev.off()
par(cex = 0.6);
par(mar = c(0,4,2,0))
par(mfrow=c(1,2))

#pick 2
plot(C9_sampleTree, main = "C9orf72", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
plot(CHMP2B_sampleTree, main = "CHMP2B", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
plot(sALS_sampleTree, main = "sALS", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
plot(FTLD_sampleTree, main = "FTLD", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
plot(VCP_sampleTree, main = "VCP", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

