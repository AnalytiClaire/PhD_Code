##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####

library (pathprint)
options(stringsAsFactors = FALSE)

thres <-225

####C9_LCM ######
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/C9orf72_LCM") #set working directory to location of data
exp_C9.LCM <- read.csv ("eset_NineP_150612_exprs.csv", header=TRUE) #assign the .csv file to a variable, column headers are true
row.names (exp_C9.LCM) <- exp_C9.LCM[,1] #specify that first column contains gene names
exp_C9.LCM<- exp_C9.LCM[,2:12] #specify that all other columns are gene expression data

C9.LCM_pathprint <- exprs2fingerprint(exp_C9.LCM, platform = "GPL570", species="human", progressBar=T) #takes the gene expression values and converts into ternary score (-1,0,1) #platform = microarray platform GEO ID

d <- apply (C9.LCM_pathprint[,1:8], 1,mean ) #d = disease, average score of each gene across all samples
c <-  apply (C9.LCM_pathprint[,9:11], 1,mean ) #c = control, average score of each gene across all samples
t <- d-c #subtract mean disease score from mean control score to find difference
C9t1<- t[order(abs(t), decreasing=T)] #order differential expression in decreasing order
c9.lcm <- (names(t1))[1:thres] #take top 'thres' values



####CHMP2B_LCM ######
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/CHMP2B")
exp_CHMP2B.LCM <- read.csv ("eset_CHMP2B_250615_exprs nooutlier.csv", header=TRUE)
row.names (exp_CHMP2B.LCM) <- exp_CHMP2B.LCM[,1]
exp_CHMP2B.LCM<- exp_CHMP2B.LCM[,2:10]

CHMP2B.LCM_pathprint <- exprs2fingerprint (exp_CHMP2B.LCM, platform = "GPL570", species="human", progressBar=T)

c <- apply (CHMP2B.LCM_pathprint[,4:9], 1,mean )
d <-  apply (CHMP2B.LCM_pathprint[,1:3], 1,mean )
t <- d-c
CHt1 <- t[order(abs(t), decreasing=T)]
CHMP2B.lcm <- (names(t1))[1:thres]

####sals_lcm###

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/FUS_SALS_LCM_CELfiles")
exp_SALS.LCM <- read.csv ("eset_SALS_LCM_260615_exprs.csv", header=TRUE)
row.names (exp_SALS.LCM) <- exp_SALS.LCM[,1]
exp_SALS.LCM<- exp_SALS.LCM[,2:11]

SALS.LCM_pathprint <- exprs2fingerprint (exp_SALS.LCM, platform = "GPL570", species="human", progressBar=T)

c <- apply (SALS.LCM_pathprint[,1:3], 1,mean )
d <-  apply (SALS.LCM_pathprint[,4:10], 1,mean )
t <- d-c
sALSt1 <- t[order(abs(t), decreasing=T)]
SALS.lcm <- (names(t1))[1:thres]

####FTLD###

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/FTD-U.brain")
FTLD <- read.csv ("FTLD_expr_tdp43.csv", header=TRUE)
row.names (FTLD) <- FTLD[,1]
FTLD <- FTLD[,2:25]

#GPL571 = Affymetrix Human Genome U113A 2.0 array
FTLD_pathprint <- exprs2fingerprint (FTLD, platform = "GPL571", species="human", progressBar=T)

c <- apply (FTLD_pathprint[,17:24], 1,mean )
d <-  apply (FTLD_pathprint[,1:16], 1,mean )
t <- d-c
FTLDt1 <- t[order(abs(t), decreasing=T)]
FTLD_FCx <- (names(t1))[1:thres]



####VCP###

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/VCP.myopathy")
VCP <- read.csv ("eset_VCP.myopathy_170715_exprs.csv", header=TRUE)
row.names (VCP) <- VCP[,1]
VCP <- VCP[,2:11]

VCP_pathprint <- exprs2fingerprint (VCP, platform = "GPL570", species="human", progressBar=T)

c <- apply (VCP_pathprint[,1:3], 1,mean )
d <-  apply (VCP_pathprint[,4:10], 1,mean )
t <- d-c
VCPt1 <- t[order(abs(t), decreasing=T)]
VCP.m <- (names(t1))[1:thres]


#Using diffPathways
C9fac <- c(1,1,1,1,1,1,1,1,0,0,0) #create vector assigning columns to disease or control
C9DP <- diffPathways(C9.LCM_pathprint, C9fac, 0.1)

CHfac <- c(1,1,1,0,0,0,0,0,0)
CHDP <- diffPathways(CHMP2B.LCM_pathprint, CHfac, 0.1)

sALSfac <- c(0,0,0,1,1,1,1,1,1,1)
sALSDP <- diffPathways(SALS.LCM_pathprint, sALSfac, 0.1)

FTLDfac <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
FTLDDP <- diffPathways(FTLD_pathprint, FTLDfac, 0.1)

VCPfac <- c(0,0,0,1,1,1,1,1,1,1)
VCPDP <- diffPathways(VCP_pathprint, VCPfac, 0.1)


#Intersect
overlap <- Reduce(intersect, list(C9DP, CHDP, sALSDP, FTLDDP, VCPDP)) #selects pathways that are present in all data sets listed
print(overlap)

setwd ("/Users/clairegreen/Desktop/")

write.csv(overlap, file = "overlap.csv")

#Heatmap

#After running pathprint
C9t1 <- as.data.frame(C9t1) 
CHt1 <- as.data.frame(CHt1)
sALSt1 <- as.data.frame(sALSt1)
FTLDt1 <- as.data.frame(FTLDt1)
VCPt1 <- as.data.frame(VCPt1)

colnames(C9t1) <- "C9Expression"
colnames(CHt1) <- "CHExpression"
colnames(sALSt1) <- "sALSExpression"
colnames(FTLDt1) <- "FTLDExpression"
colnames(VCPt1) <- "VCPExpression"

C9t1[,2] <- rownames(C9t1)
CHt1[,2] <- rownames(CHt1)
sALSt1[,2] <- rownames(sALSt1)
FTLDt1[,2] <- rownames(FTLDt1)
VCPt1[,2] <- rownames(VCPt1)

all.fingerprints <- merge(x = C9t1, y = CHt1, by.x= "V2", by.y="V2")
all.fingerprints <- merge(x = all.fingerprints, y = sALSt1, by.x= "V2", by.y="V2")
all.fingerprints <- merge(x = all.fingerprints, y = FTLDt1, by.x= "V2", by.y="V2")
all.fingerprints <- merge(x = all.fingerprints, y = VCPt1, by.x= "V2", by.y="V2")


all <- C9.LCM_pathprint


library(data.table)
setDT(all, keep.rownames = TRUE)[]

overlap <- as.data.frame(overlap) #conver to data frame

sig.names <- read.delim(file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/DEPathways.txt")
sig.fingerprints <- merge(overlap, C9all, by.x = "overlap", by.y = "rn" ) #take only pathways implicated by intersection

rownames(sig.fingerprints) <- sig.fingerprints[,1]
sig.fingerprints[,1] <- NULL

sig.fingerprints <- as.matrix(sig.fingerprints) #must be converted to numeric data frame for heatmap to work

heatmap(sig.fingerprints)
