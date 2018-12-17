
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression")
C9 <- read.csv("C9rankeduniqueresult.csv")
sals <- read.csv("sALSrankeduniqueresult.csv")
ftld <- read.csv("FTLDrankeduniqueresult.csv")
GRNftld <- read.csv("FTLD_GRNrankeduniqueresult.csv")
sftld <- read.csv("FTLD_SPrankeduniqueresult.csv")
vcp <- read.csv("VCPrankeduniqueresult.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_results_keepfiltering.csv")
pet <- pet[!duplicated(pet$hgnc_symbol),]
petC9 <- read.csv("PET_C9_results_keepfiltering.csv")
petC9 <- petC9[!duplicated(pet$hgnc_symbol),]
petsALS <- read.csv("PET_sALS_results_keepfiltering.csv")
petsALS <- petsALS[!duplicated(pet$hgnc_symbol),]

rav <- read.csv("RAV_results_keepfiltering.csv")
rav <- rav[!duplicated(rav$hgnc_symbol),]


genelist <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NuggetGenes.txt")


subsetC9 <- subset(C9, C9$Gene.Symbol %in% genelist, drop = TRUE)
subsetC9 <- subsetC9[!duplicated(subsetC9$Gene.Symbol),]
rownames(subsetC9) <- subsetC9$Gene.Symbol
subsetC9 <- subsetC9[order(row.names(subsetC9)),]

subsetsals <- subset(sals, sals$Gene.Symbol %in% genelist, drop = TRUE)
subsetsals <- subsetsals[!duplicated(subsetsals$Gene.Symbol),]
rownames(subsetsals) <- subsetsals$Gene.Symbol
subsetsals <- subsetsals[order(row.names(subsetsals)),]

subsetftld <- subset(ftld, ftld$Gene.Symbol %in% genelist, drop = TRUE)
subsetftld <- subsetftld[!duplicated(subsetftld$Gene.Symbol),]
rownames(subsetftld) <- subsetftld$Gene.Symbol
subsetftld <- subsetftld[order(row.names(subsetftld)),]

subsetsftld <- subset(sftld, sftld$Gene.Symbol %in% genelist, drop = TRUE)
subsetsftld <- subsetsftld[!duplicated(subsetsftld$Gene.Symbol),]
rownames(subsetsftld) <- subsetsftld$Gene.Symbol
subsetsftld <- subsetsftld[order(row.names(subsetsftld)),]

subsetGRNftld <- subset(GRNfltd, GRNftld$Gene.Symbol %in% genelist, drop = TRUE)
subsetGRNftld <- subsetGRNftld[!duplicated(subsetGRNftld$Gene.Symbol),]
rownames(subsetGRNftld) <- subsetGRNftld$Gene.Symbol
subsetGRNftld <- subsetGRNftld[order(row.names(subsetGRNftld)),]

subsetvcp <- subset(vcp, vcp$Gene.Symbol %in% genelist, drop = TRUE)
subsetvcp <- subsetvcp[!duplicated(subsetvcp$Gene.Symbol),]
rownames(subsetvcp) <- subsetvcp$Gene.Symbol
subsetvcp <- subsetvcp[order(row.names(subsetvcp)),]

subsetpet <- subset(pet, pet$hgnc_symbol %in% genelist, drop = TRUE)
subsetpet <- subsetpet[!duplicated(subsetpet$hgnc_symbol),]
rownames(subsetpet) <- subsetpet$Gene.Symbol
subsetpet <- subsetpet[order(subsetpet$hgnc_symbol),]
rownames(subsetpet) <- subsetpet$hgnc_symbol

subsetpetC9 <- subset(petC9, petC9$hgnc_symbol %in% genelist, drop = TRUE)
subsetpetC9 <- subsetpetC9[!duplicated(subsetpetC9$hgnc_symbol),]
rownames(subsetpetC9) <- subsetpetC9$Gene.Symbol
subsetpetC9 <- subsetpetC9[order(subsetpetC9$hgnc_symbol),]
rownames(subsetpetC9) <- subsetpetC9$hgnc_symbol

subsetpetsALS <- subset(petsALS, petsALS$hgnc_symbol %in% genelist, drop = TRUE)
subsetpetsALS <- subsetpetsALS[!duplicated(subsetpetsALS$hgnc_symbol),]
rownames(subsetpetsALS) <- subsetpetsALS$Gene.Symbol
subsetpetsALS <- subsetpetsALS[order(subsetpetsALS$hgnc_symbol),]
rownames(subsetpetsALS) <- subsetpetsALS$hgnc_symbol

subsetrav <- subset(rav, rav$hgnc_symbol %in% genelist, drop = TRUE)
subsetrav <- subsetrav[!duplicated(subsetrav$hgnc_symbol),]
rownames(subsetrav) <- subsetrav$Gene.Symbol
subsetrav <- subsetrav[order(subsetrav$hgnc_symbol),]
rownames(subsetrav) <- subsetrav$hgnc_symbol


#PatientvsControl Summary Data
C9gen <- subsetC9[,14:20]
salsgen <- subsetsals[,14:20]
ftldgen <- subsetftld[,13:19]
sfltdgen <- subsetsftld[,13:19]
GRNfltdgen <- subsetGRNftld[,13:19]
vcpgen <- subsetvcp[,14:20]
petgen <- subsetpet[,c(4:9,37)]
petC9gen <- subsetpet[,c(4:9,37)]
petsalsgen <- subsetpet[,c(4:9,37)]
ravgen <- subsetrav[,c(4:9,32)]

##Patients and control Expression data
# C9gen <- subsetC9[,24:31]
# salsgen <- subsetsals[,24:30]
# ftldgen <- subsetftld[,28:43]
# sfltd <- subsetftld[,34:43]
# GRNfltd <- subsetftld[,28:33]
# vcpgen <- subsetvcp[,24:30]
# petgen <- subsetpet[19:35]
# petC9 <- subsetpet[19:25]
# petsals <- subsetpet[26:35]
# ravgen <- subsetrav[,18:30]

##Patients only Expression data
# C9gen <- subsetC9[,21:31]
# salsgen <- subsetsals[,21:30]
# ftldgen <- subsetftld[,20:43]
# sfltd <- subsetftld[,c(20:27,34:43)]
# GRNfltd <- subsetftld[,c(20:27,28:33)]
# vcpgen <- subsetvcp[,21:30]
# petgen <- subsetpet[,10:35]
# petC9 <- subsetpet[,10:25]
# petsals <- subsetpet[,c(10:18,26:35)]
# ravgen <- subsetrav[,10:30]

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/IPA/")
write.csv(C9gen, "C9_nugget.csv")
write.csv(salsgen, "sals_nugget.csv")
write.csv(ftldgen, "ftld_nugget.csv")
write.csv(sfltdgen, "sftld_nugget.csv")
write.csv(GRNfltdgen, "GRNftld_nugget.csv")
write.csv(vcpgen, "vcp_nugget.csv")
write.csv(petgen, "pet_nugget.csv")
write.csv(petC9gen, "petC9_nugget.csv")
write.csv(petsalsgen, "petsals_nugget.csv")
write.csv(ravgen, "rav_nugget.csv")



data_all <- cbind(C9gen, salsgen, ftldgen, vcpgen, petgen, ravgen)

write.csv(LFC, "LogFoldChange_Nugget_perDataset.csv", row.names = F, quote = F)
