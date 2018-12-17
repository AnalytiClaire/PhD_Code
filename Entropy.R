library(entropy)

#Run TDP43_signature first

x=getEntropy(C9.LCM_pathprint,1)
x.selectC9=C9.LCM_pathprint[names(x)[which(x>0.5)],]
heatmap(x.selectC9)

x=getEntropy(CHMP2B.LCM_pathprint,1)
x.selectCH=CHMP2B.LCM_pathprint[names(x)[which(x>0.5)],]
heatmap(x.selectCH)

x=getEntropy(SALS.LCM_pathprint,1)
x.selectsals=SALS.LCM_pathprint[names(x)[which(x>0.5)],]
heatmap(x.selectsals)

x=getEntropy(FTLD_pathprint,1)
x.selectFTLD=FTLD_pathprint[names(x)[which(x>0.5)],]
heatmap(x.selectFTLD)

x=getEntropy(VCP_pathprint,1)
x.selectVCP=VCP_pathprint[names(x)[which(x>0.5)],]
heatmap(x.selectVCP)

#Using diffPathways

DEthresh <- 0

C9fac <- c(1,1,1,1,1,1,1,1,0,0,0) #create vector assigning columns to disease or control
C9DP <- diffPathways(x.selectC9, C9fac, DEthresh)

CHfac <- c(1,1,1,0,0,0,0,0,0)
CHDP <- diffPathways(x.selectCH, CHfac, DEthresh)

sALSfac <- c(0,0,0,1,1,1,1,1,1,1)
sALSDP <- diffPathways(x.selectsals, sALSfac, DEthresh)

FTLDfac <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
FTLDDP <- diffPathways(x.selectFTLD, FTLDfac, DEthresh)

VCPfac <- c(0,0,0,1,1,1,1,1,1,1)
VCPDP <- diffPathways(x.selectVCP, VCPfac, DEthresh)


#Intersect
overlap <- Reduce(intersect, list(C9DP, CHDP, sALSDP, FTLDDP, VCPDP)) #selects pathways that are present in all data sets listed
print(overlap)
