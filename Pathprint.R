library (pathprint)
options(stringsAsFactors = FALSE)

####C9_LCM ######
setwd ("/Users/clairegreen/Documents/PhD/Data Sets/C9orf72_LCM")
exp_C9.LCM <- read.csv ("eset_NineP_150612_exprs.csv", header=TRUE)
row.names (exp_C9.LCM) <- exp_C9.LCM[,1]
exp_C9.LCM<- exp_C9.LCM[,2:12]

C9.LCM_pathprint <- exprs2fingerprint (exp_C9.LCM, platform = "GPL570", species="human", progressBar=T)

d <- apply (C9.LCM_pathprint[,1:8], 1,mean )
c <-  apply (C9.LCM_pathprint[,9:11], 1,mean )
t <- d-c
t1 <- t[order(abs(t), decreasing=T)]
c9.lcm <- (names(t1))[1:100]




# ####CHMP2B_LCM ######
# setwd ("/Users/clairegreen/Documents/PhD/Data Sets/CHMP2B")
# exp_CHMP2B.LCM <- read.csv ("eset_CHMP2B_250615_exprs.csv", header=TRUE)
# row.names (exp_CHMP2B.LCM) <- exp_CHMP2B.LCM[,1]
# exp_CHMP2B.LCM<- exp_CHMP2B.LCM[,2:11]
# 
# CHMP2B.LCM_pathprint <- exprs2fingerprint (exp_CHMP2B.LCM, platform = "GPL570", species="human", progressBar=T)
# 
# c <- apply (CHMP2B.LCM_pathprint[,1:7], 1,mean )
# d <-  apply (CHMP2B.LCM_pathprint[,8:10], 1,mean )
# t <- d-c
# t1 <- t[order(abs(t), decreasing=T)]
# CHMP2B.lcm <- (names(t1))[1:100]





####sals_lcm###

setwd ("/Users/clairegreen/Documents/PhD/Data Sets/FUS_SALS_LCM_CELfiles")
exp_SALS.LCM <- read.csv ("eset_SALS_LCM_260615_exprs.csv", header=TRUE)
row.names (exp_SALS.LCM) <- exp_SALS.LCM[,1]
exp_SALS.LCM<- exp_SALS.LCM[,2:11]

SALS.LCM_pathprint <- exprs2fingerprint (exp_SALS.LCM, platform = "GPL570", species="human", progressBar=T)

c <- apply (SALS.LCM_pathprint[,1:3], 1,mean )
d <-  apply (SALS.LCM_pathprint[,4:10], 1,mean )
t <- d-c
t1 <- t[order(abs(t), decreasing=T)]
SALS.lcm <- (names(t1))[1:100]

####FTLD###

setwd ("/Users/clairegreen/Documents/PhD/Data Sets/FTD-U.brain")
FTLD <- read.csv ("eset_FTD.U.brain_170715_exprs.csv",header=TRUE)
row.names (FTLD) <- FTLD[,1]
FTLD<- FTLD[,2:57]

FCx.Con <- c(1,4,5,7,8,11,13,15)
FCx.PRGN <- c(18,20,22,24,27,30)
FCx.SFTD <- c(33,35,38,41,44,45,48,50,52,55)

FTLD_pathprint_<- exprs2fingerprint (FTLD, platform = "GPL570",species="human",progressBar=T)

c <- apply (FTLD_pathprint[,c(FCx.Con)], 1,mean )
d <-  apply (FTLD_pathprint[,c(FCx.PRGN, FCx.SFTD)], 1,mean )
t <- d-c
t1 <- t[order(abs(t), decreasing=T)]
FTLD_FCx <- (names(t1))[1:100]

####VCP###

setwd ("/Users/clairegreen/Documents/PhD/Data Sets/VCP.myopathy")
VCP <- read.csv ("eset_VCP.myopathy_170715_exprs.csv", header=TRUE)
row.names (VCP) <- VCP[,1]
VCP<- VCP[,2:11]

VCP_pathprint <- exprs2fingerprint (VCP, platform = "GPL570", species="human", progressBar=T)

c <- apply (VCP_pathprint[,1:3], 1,mean )
d <-  apply (VCP_pathprint[,4:10], 1,mean )
t <- d-c
t1 <- t[order(abs(t), decreasing=T)]
VCP <- (names(t1))[1:100]

#### intersect

i1 <- intersect (c9.lcm, CHMP2B.lcm)
print (i1)

i2 <- intersect (CHMP2B.lcm, SALS.lcm)
print (i2)

overlap_lcm <- intersect (i1, i2)

i3 <- intersect (FTLD_FCx, overlap_lcm)

i4 <- intersect (i3,VCP)

print (overlap_lcm)


