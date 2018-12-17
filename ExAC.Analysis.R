#enrichment permutation#

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/ExAC/")
#Load ExAC Data
Exac.All <- read.table(file = "fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", header = TRUE)
Exac.95 <- read.table(file = "exac.pli.0.95.txt", header = TRUE)
Exac <- Exac.All[c("gene", "pLI", "mis_z")]
rownames(Exac) <- Exac$gene
Exac[,1] <- NULL

#Extract my genes from that list
DEG.ExAC <- Exac[c("ACTN1", "ANXA1", "BBIP1", "BGN", "BPTF", "CDH11", "CREG1", "CSRP1", "CST3", "DCN", "GBAS", "JAG1",
                   "KCTD12", "KPNA6", "MPHOSPH9", "MXI1", "NDUFS5 /// RPL10", "NKTR", "NUTF2 /// NUTF2P4", "LOC101927673 /// OTUB1", 
                   "PFDN1", "PLEKHB1", "PLOD2", "POGZ", "PPP1R7", "PPP2CA", "PPP2CB", "PRPF3", "PTPN13", "RAB40B", "RPL35A", "RPL37", 
                   "SCN1B", "SERBP1", "SF3B1", "SPARC", "STMN1", "TARDBP", "TCF4", "TUBB4B", "TUG1", "VPS13B", "ZFP36", "ZFYVE26", "ZNF518A"),]

#Extract top 5% pLI from that list

Top5 <- Exac.All[Exac.All$gene %in% l,]
DEG.Top5 <- Top5[Top5$gene %in% g,]
DEG45 <- Top5[Top5$gene %in% o1,]

DEGtotal.pLI <- sum(DEG45$pLI)
DEGtotal.mis_Z <- sum(DEG45$mis_z)


test <- DEGtotal.pLI
m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

for (j in 1:m)
{
  random <- sample (Top5$pLI, size=26, replace=FALSE)
  r[j] <- sum(random)
}

test1 <- which(r > test)  # count number of times r is larger than test value
result <- (length(test1)/m) # calculate P value
mean(r)




test <- DEGtotal.mis_Z
m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

for (j in 1:m)
{
  random <- sample (Top5$mis_z, size=13, replace=FALSE)
  r[j] <- sum(random)
}

test1 <- which(r > test)  # count number of times r is larger than test value
result <- (length(test1)/m) # calculate P value
mean(r)






DEG.Top5 <- Top5[c("ACTN1" , "BPTF"  ,"CDH11"  ,"JAG1"  , "KPNA6" , "NKTR"   ,"POGZ"   ,"PPP2CA" ,"PRPF3"  ,"SERBP1" ,"SF3B1"  ,"TARDBP" ,"TCF4"),]



genes <- read.table(file = "ExAC5.DEGs.txt")
genes <- genes$V1

library(biomaRt)
genes <- rownames(DEG.Top5)
mart <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset="hsapiens_snp")

mart_back <- getBM(attributes =c("refsnp_id"), filters="feature_stable_id", values=genes,  mart=mart)
  