diffPathwayTable <- function (fingerprints, fac, threshold) 
{
  fac <- as.factor(fac)
  levels <- levels(fac)
  if (length(levels) != 2) 
    stop("fac must contain 2 levels")
  if (length(fac) != ncol(fingerprints)) 
    stop("fac length must equal nubmer of fingerprints")
  if ((threshold < 0 | threshold > 2)) 
    stop("threshold value should be between 0 and 2")
  if (sum(fac == levels[1]) == 1) {
    print("N.B. only 1 sample for group 1")
    mean.1 <- fingerprints[, fac == levels[1]]
  }
  else {
    mean.1 <- apply(fingerprints[, fac == levels[1]], 1, 
                    mean)
  }
  if (sum(fac == levels[2]) == 1) {
    print("N.B. only 1 sample for group 2")
    mean.2 <- fingerprints[, fac == levels[2]]
  }
  else {
    mean.2 <- apply(fingerprints[, fac == levels[2]], 1, 
                    mean)
  }
  mean.diff <- mean.2 - mean.1
  meandiff <- abs(mean.diff)
  signifPathawys <- na.omit(names(mean.diff)[abs(mean.diff) > 
                                               threshold])
  return(as.data.frame(mean.diff))
  return(as.vector(signifPathawys))
}

##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####

##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####

library (pathprint)
options(stringsAsFactors = FALSE)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/NormalisedExpressionMatrices/") #set working directory to location of data

####C9_LCM ######
exp_C9.LCM <- read.csv ("C9eset.csv", header=TRUE, row.names = 1) #assign the .csv file to a variable, column headers are true
C9.LCM_pathprint <- exprs2fingerprint(exp_C9.LCM, platform = "GPL570", species="human", progressBar=T) #takes the gene expression values and converts into ternary score (-1,0,1) #platform = microarray platform GEO ID
vec.c9 <- c(0,0,0,1,1,1,1,1,1,1,1)


# ####CHMP2B_LCM ######
# # setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/CHMP2B")
# exp_CHMP2B.LCM <- read.csv ("CHeset.csv", header=TRUE, row.names = 1)
# # row.names (exp_CHMP2B.LCM) <- exp_CHMP2B.LCM[,1]
# # exp_CHMP2B.LCM<- exp_CHMP2B.LCM[,2:10]
# 
# CHMP2B.LCM_pathprint <- exprs2fingerprint (exp_CHMP2B.LCM, platform = "GPL570", species="human", progressBar=T)
# vec.ch <- c(0,0,0,0,0,0,1,1,1)


####sals_lcm###

exp_SALS.LCM <- read.csv ("sALSeset.csv", header=TRUE, row.names = 1)
SALS.LCM_pathprint <- exprs2fingerprint (exp_SALS.LCM, platform = "GPL570", species="human", progressBar=T)
vec.sals <- c(0,0,0,1,1,1,1,1,1,1)


####FTLD###
exp_FTLD <- read.csv ("FTLDeset.csv", header=TRUE, row.names = 1)
FTLD_pathprint <- exprs2fingerprint (FTLD, platform = "GPL571", species="human", progressBar=T)
vec.FTLD <- c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)


####VCP###
exp_VCP <- read.csv ("VCPeset.csv", header=TRUE, row.names = 1)
VCP_pathprint <- exprs2fingerprint (VCP, platform = "GPL570", species="human", progressBar=T)
vec.vcp <- c(0,0,0,1,1,1,1,1,1,1)



c9.lcm <- diffPathwayTable(C9.LCM_pathprint, vec.c9, thres)
# CHMP2B.lcm <- diffPathwayTable(CHMP2B.LCM_pathprint, vec.ch, thres)
SALS.lcm <- diffPathwayTable(SALS.LCM_pathprint, vec.sals, thres)
FTLD_FCx <- diffPathwayTable(FTLD_pathprint, vec.FTLD, thres)
VCP.m <- diffPathwayTable(VCP_pathprint, vec.vcp, thres)

# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/PathprintThreshold")
# write.csv(c9.lcm, "C9.csv")
# write.csv(SALS.lcm, "sALS.csv")
# write.csv(CHMP2B.lcm, "CHMP2B.csv")
# write.csv(FTLD_FCx, "FTLD.csv")
# write.csv(VCP.m, "VCP.csv")


library(data.table)
setDT(c9.lcm, keep.rownames = TRUE)
# setDT(CHMP2B.lcm, keep.rownames = TRUE)
setDT(SALS.lcm, keep.rownames = TRUE)
setDT(FTLD_FCx, keep.rownames = TRUE)
setDT(VCP.m, keep.rownames = TRUE)

c9.lcm <- c9.lcm[order(c9.lcm[,2], decreasing = TRUE),]
# CHMP2B.lcm <- CHMP2B.lcm[order(CHMP2B.lcm$meandiff, decreasing = TRUE),]
SALS.lcm <- SALS.lcm[order(SALS.lcm$mean.diff, decreasing = TRUE),]
FTLD_FCx <- FTLD_FCx[order(FTLD_FCx$mean.diff, decreasing = TRUE),]
VCP.m <- VCP.m[order(VCP.m$mean.diff, decreasing = TRUE),]

# thresh <- 200
# 
# c9thresh <- c9.lcm[1:thresh,]
# chthresh <- CHMP2B.lcm[1:thresh,]
# salsthresh <- SALS.lcm[1:thresh,]
# ftldthresh <- FTLD_FCx[1:thresh,]
# vcpthresh <- VCP.m[1:thresh,]



overlap <- Reduce(intersect, list(c9thresh$rn, chthresh$rn, salsthresh$rn, ftldthresh$rn, vcpthresh$rn)) #selects pathways that are present in all data sets listed
overlap

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/PathprintThreshold")
write.table(overlap, "Pathprint_thresholdrevision.txt", quote = F, row.names = F, col.names = F)



pathmatrix <- data.frame(path=c9.lcm$rn, 
                    C9orf72=c9.lcm$mean.diff,
                    sALS=SALS.lcm$mean.diff, 
                    FTLD=FTLD_FCx$mean.diff,
                    VCP=VCP.m$mean.diff)


######
rownames(pathmatrix) <- pathmatrix$path
pathmatrix[,1] <- NULL

upthres <- 0.5
pathup <- subset(pathmatrix, C9orf72 > upthres & sALS > upthres & FTLD > upthres & VCP > upthres)

downthres <- -0.5
pathdown <- subset(pathmatrix, C9orf72 < downthres & sALS < downthres & FTLD < downthres & VCP < downthres)

setwd("")
write.csv(pathup, "uppathways_0.5.csv")
