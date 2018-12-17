setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/")

C9cor <- read.csv("C9orf72_coexpression.csv")
C9cor <- C9cor[order(C9cor$X),]
salscor <- read.csv("sALS_coexpression.csv")
salscor <- salscor[order(salscor$X),]
FTLDcor <- read.csv("FTLD_coexpression.csv")
FTLDcor <- FTLDcor[order(FTLDcor$X),]
VCPcor <- read.csv("VCP_coexpression.csv")
VCPcor <- VCPcor[order(VCPcor$X),]
PETcor <- read.csv("PET_coexpression.csv")
PETcor <- PETcor[order(PETcor$X),]
RAVcor <- read.csv("RAV_coexpression.csv")
RAVcor <- RAVcor[order(RAVcor$X),]

cor <- data.frame(row.names = C9cor$X,
                     C9orf72 = C9cor$reg.mat,
                     sALS = salscor$reg.mat,
                     FTLD = FTLDcor$reg.mat,
                     VCP = VCPcor$reg.mat,
                     PET = PETcor$reg.mat,
                     RAV = RAVcor$reg.mat)

poscor <- cor[apply(cor, MARGIN = 1, function(x) all(x > 0)),]
negcor <- cor[apply(cor, MARGIN = 1, function(x) all(x < 0)),]
all_cor <- rbind(poscor, negcor)

all_cor$mean <- apply(all_cor, 1, mean)

all_cor_mean <- as.data.frame(all_cor$mean, row.names = rownames(all_cor))
all_cor_mean$INT <- "Correlation"
write.csv(all_cor_mean, "MeanCorrelation_TDP43.csv")
