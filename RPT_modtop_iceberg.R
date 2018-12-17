#### Random Permutation Test for Edge Connectivity ####

library(gdata)
library(tictoc)


C9 <- read.csv("C9uniquegene_samples.csv", row.names = 1)
C9 <- C9[,4:11]
VCP <- read.csv("VCPuniquegene_samples.csv", row.names = 1)
VCP <- VCP[,4:10]
FTLD <- read.csv("FTLDuniquegene_samples.csv", row.names = 1)
FTLD <- FTLD[,9:24]
sALS <- read.csv("sALSuniquegene_samples.csv", row.names = 1)
sALS <- sALS[,4:10]
PET <- read.csv("PET_results_keepfiltering.csv")
rownames(PET) <- PET$hgnc_symbol
PET <- PET[,19:35]
RAV <- read.csv("RAV_results_keepfiltering.csv")
rownames(RAV) <- RAV$hgnc_symbol
RAV <- RAV[,18:30]

#Find the gene names that all datasets have in common
comgene <- Reduce(intersect, list(rownames(C9), rownames(VCP), rownames(FTLD),
                                  rownames(sALS), rownames(PET), rownames(RAV)))

m=10 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

tic()
for (k in 1:m){
  #pick 72 random genes (for dummy module)
  modsize = 72
  RanModule <- sample(comgene, size = modsize, replace = FALSE)
  
  #Extract PPI network genes from each dataset
  
  C9.72 <- subset(C9, rownames(C9) %in% RanModule)
  VCP.72 <- subset(VCP, rownames(VCP) %in% RanModule)
  FTLD.72 <- subset(FTLD, rownames(FTLD) %in% RanModule)
  sALS.72 <- subset(sALS, rownames(sALS) %in% RanModule)
  PET.72 <- subset(PET, rownames(PET) %in% RanModule)
  RAV.72 <- subset(RAV, rownames(RAV) %in% RanModule)
  
  ###### C9 #########
  CorExprMat <- t(C9.72)
  test <- CorExprMat
  reg <- matrix(0, ncol(test), ncol(test))
  for (i in 1:ncol(test)){
    for (j in 1:ncol(test)){
      reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
    }}
  rownames(reg) <- colnames(reg) <- colnames(test)
  ptri <- reg
  ptri[lower.tri(reg, diag = TRUE)] <- NA
  p.vec <- unmatrix(ptri)
  p.vec <- na.omit(p.vec)
  p.mat.C9 <- as.data.frame(p.vec)
  colnames(p.mat.C9) <- "C9"
  p.mat.C9$names <- rownames(p.mat.C9)
  p.mat.C9 <- p.mat.C9[order(p.mat.C9$names),]
  
  
  ###### sALS #########
  CorExprMat <- t(sALS.72)
  test <- CorExprMat
  reg <- matrix(0, ncol(test), ncol(test))
  for (i in 1:ncol(test)){
    for (j in 1:ncol(test)){
      reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
    }}
  rownames(reg) <- colnames(reg) <- colnames(test)
  ptri <- reg
  ptri[lower.tri(reg, diag = TRUE)] <- NA
  p.vec <- unmatrix(ptri)
  p.vec <- na.omit(p.vec)
  p.mat.sALS <- as.data.frame(p.vec)
  colnames(p.mat.sALS) <- "sALS"
  p.mat.sALS$names <- rownames(p.mat.sALS)
  p.mat.sALS <- p.mat.sALS[order(p.mat.sALS$names),]
  
  ###### FTLD #########
  CorExprMat <- t(FTLD.72)
  test <- CorExprMat
  reg <- matrix(0, ncol(test), ncol(test))
  for (i in 1:ncol(test)){
    for (j in 1:ncol(test)){
      reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
    }}
  rownames(reg) <- colnames(reg) <- colnames(test)
  ptri <- reg
  ptri[lower.tri(reg, diag = TRUE)] <- NA
  p.vec <- unmatrix(ptri)
  p.vec <- na.omit(p.vec)
  p.mat.FTLD <- as.data.frame(p.vec)
  colnames(p.mat.FTLD) <- "FTLD"
  p.mat.FTLD$names <- rownames(p.mat.FTLD)
  p.mat.FTLD <- p.mat.FTLD[order(p.mat.FTLD$names),]
  
  ###### VCP #########
  CorExprMat <- t(VCP.72)
  test <- CorExprMat
  reg <- matrix(0, ncol(test), ncol(test))
  for (i in 1:ncol(test)){
    for (j in 1:ncol(test)){
      reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
    }}
  rownames(reg) <- colnames(reg) <- colnames(test)
  ptri <- reg
  ptri[lower.tri(reg, diag = TRUE)] <- NA
  p.vec <- unmatrix(ptri)
  p.vec <- na.omit(p.vec)
  p.mat.VCP <- as.data.frame(p.vec)
  colnames(p.mat.VCP) <- "VCP"
  p.mat.VCP$names <- rownames(p.mat.VCP)
  p.mat.VCP <- p.mat.VCP[order(p.mat.VCP$names),]
  
  ###### PET #########
  CorExprMat <- t(PET.72)
  test <- CorExprMat
  reg <- matrix(0, ncol(test), ncol(test))
  for (i in 1:ncol(test)){
    for (j in 1:ncol(test)){
      reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
    }}
  rownames(reg) <- colnames(reg) <- colnames(test)
  ptri <- reg
  ptri[lower.tri(reg, diag = TRUE)] <- NA
  p.vec <- unmatrix(ptri)
  p.vec <- na.omit(p.vec)
  p.mat.PET <- as.data.frame(p.vec)
  colnames(p.mat.PET) <- "PET"
  p.mat.PET$names <- rownames(p.mat.PET)
  p.mat.PET <- p.mat.PET[order(p.mat.PET$names),]
  
  ###### RAV #########
  CorExprMat <- t(RAV.72)
  test <- CorExprMat
  reg <- matrix(0, ncol(test), ncol(test))
  for (i in 1:ncol(test)){
    for (j in 1:ncol(test)){
      reg[i,j] <- cor.test(test[,i], test[,j], method = "spearman")$estimate
    }}
  rownames(reg) <- colnames(reg) <- colnames(test)
  ptri <- reg
  ptri[lower.tri(reg, diag = TRUE)] <- NA
  p.vec <- unmatrix(ptri)
  p.vec <- na.omit(p.vec)
  p.mat.RAV <- as.data.frame(p.vec)
  colnames(p.mat.RAV) <- "RAV"
  p.mat.RAV$names <- rownames(p.mat.RAV)
  p.mat.RAV <- p.mat.RAV[order(p.mat.RAV$names),]
  
  C9_cor.5 <- p.mat.C9[p.mat.C9$C9 > 0.5 | p.mat.C9$C9 < -0.5,]
  sALS_cor.5 <- p.mat.sALS[p.mat.sALS$sALS > 0.5 | p.mat.sALS$sALS < -0.5,]
  FTLD_cor.5 <- p.mat.FTLD[p.mat.FTLD$FTLD > 0.5 | p.mat.FTLD$FTLD < -0.5,]
  VCP_cor.5 <- p.mat.VCP[p.mat.VCP$VCP > 0.5 | p.mat.VCP$VCP < -0.5,]
  PET_cor.5 <- p.mat.PET[p.mat.PET$PET > 0.5 | p.mat.PET$PET < -0.5,]
  RAV_cor.5 <- p.mat.RAV[p.mat.RAV$RAV > 0.5 | p.mat.RAV$RAV < -0.5,]
  
  
  Commonedge <- Reduce(intersect, list(C9_cor.5$names, sALS_cor.5$names, FTLD_cor.5$names,
                                       VCP_cor.5$names, PET_cor.5$names, RAV_cor.5$names))
  
  
  #Subset each dataset with these common names so they are all the same size
  C9_CE <- subset(C9_cor.5, C9_cor.5$names %in% Commonedge)
  C9_CE <- C9_CE[order(C9_CE$names),]
  sALS_CE <- subset(sALS_cor.5, sALS_cor.5$names %in% Commonedge)
  sALS_CE <- sALS_CE[order(sALS_CE$names),]
  FTLD_CE <- subset(FTLD_cor.5, FTLD_cor.5$names %in% Commonedge)
  FTLD_CE <- FTLD_CE[order(FTLD_CE$names),]
  VCP_CE <- subset(VCP_cor.5, VCP_cor.5$names %in% Commonedge)
  VCP_CE <- VCP_CE[order(VCP_CE$names),]
  PET_CE <- subset(PET_cor.5, PET_cor.5$names %in% Commonedge)
  PET_CE <- PET_CE[order(PET_CE$names),]
  RAV_CE <- subset(RAV_cor.5, RAV_cor.5$names %in% Commonedge)
  RAV_CE <- RAV_CE[order(RAV_CE$names),]
  
  
  CommonGroup <- data.frame(row.names = C9_CE$names,
                            C9 = C9_CE$C9,
                            sALS = sALS_CE$sALS,
                            FTLD = FTLD_CE$FTLD,
                            VCP = VCP_CE$VCP,
                            PET = PET_CE$PET,
                            RAV = RAV_CE$RAV)
  
  
  CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
  CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]
  
  CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
  CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
  
  r[k] <- length(rownames(CG_samedir))
}
toc()
r
test1 <- which(r > test)  # count number of times r is larger than test value
result <- sum((length(test1)+1))/(m+1) # calculate P value
result
mean <- mean(r)
range <- paste(r[1], r[2], sep = ":")

table <- data.frame(NumOverTest = length(test1),
                    Pval = result,
                    mean = mean,
                    range = range)

write.csv(table, "edgeenrichmentresults.csv", row.names = F)
