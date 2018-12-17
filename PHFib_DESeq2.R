#### Analysis of PH's TDP-43 Fibroblast data ####
options(scipen=1)

library(DESeq2)
library(biomaRt)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/PH_Fibroblasts/")
counts <- read.table("Counts/annotated_combinedcounts.txt", header = T, row.names = 1)


counts.g<- counts[,1:30] #remove last column
counts.g <- counts.g[,c(1:3,7,4:6,8:13,17,14:16,18:23,27,24:26,28:30)] #rearrange columns into control/patient order


##### QC #####

# PCA 2D #
library(ggplot2)
library(gridExtra)
names <- colnames(counts.g)
names <- sapply(strsplit(names,split= "_[aA-zZ]"),'[',1)


dev.off() #removes previous plot
pca_2D <- prcomp(t(counts.g))
PC <- c(1,2)
plot(pca_2D$x[,PC], pch=18, cex=2 , col=c(rep("chartreuse",4), 
                                             rep("chartreuse4",6), 
                                             rep("dodgerblue",4),
                                             rep("dodgerblue4",6), 
                                             rep("gold2",4),
                                             rep("peru",6)),
     main = "2D PCA Plot: PC1 vs PC2") #plot disease

#add legend
legend(-300000, -210000,pch=18, 
       legend=c("Cytoplasm - Control", 
                "Cytoplasm - Patient",
                "Nucleus - Control",
                "Nucleus - Patient", 
                "WCT - Control",
                "WCT - Patient"),
       col=c("chartreuse","chartreuse4","dodgerblue","dodgerblue4","gold2","peru"), 
       cex=0.9,
       pt.cex = 2) 
                                                                 

text(pca_2D$x[,PC], labels = names,adj = 1.3, cex = 0.8, pos = 3) #label data points



#3d PCA
library(rgl)
library(pca3d)

pca <- prcomp(t(counts.g)) #conduct PCA
pca3d(pca$x[,1:3],
      col=c(rep("chartreuse",4),
            rep("chartreuse4",6), 
            rep("dodgerblue",4),
            rep("dodgerblue4",6), 
            rep("gold2",4),
            rep("peru",6)),
      new= TRUE, 
      show.plane = FALSE, 
      legend=1) #plot 3D
text3d(pca$x[,1:3], text=names, adj=1.3, color="black", cex = 0.7) #add labels

#QQ plot
qqplot(x = counts.g[,1], y = counts.g[,5], xlab = "Patient_201", ylab = "Patient_78")


## Box plots ##

dev.off()
logdata <- log(counts.g)


boxplot(counts.g, 
        las = 2,
        names = names,
        par(mar = c(7, 5, 4, 2)+ 0.1),
        ylim=c(0.1,500000),
        log = "y",
        outcol=c(rep("chartreuse",4),
                 rep("chartreuse4",6), 
                 rep("dodgerblue",4),
                 rep("dodgerblue4",6), 
                 rep("gold2",4),
                 rep("peru",6)),
        ylab = "Log Count")



##### DIFFERENTIAL EXPRESSION #####

#Select out comparisons to be made
cyt <- counts.g[,1:10]
nuc <- counts.g[,11:20]
WCT <- counts.g[,21:30]

cyt.v.nuc.pat <- data.frame(c(cyt[5:10], nuc[5:10]))
rownames(cyt.v.nuc.pat) <- rownames(cyt)
cyt.v.nuc.con <- data.frame(c(cyt[1:4], nuc[1:4]))
rownames(cyt.v.nuc.con) <- rownames(cyt)

WCT.v.nuc.pat <- data.frame(c(WCT[5:10], nuc[5:10]))
rownames(WCT.v.nuc.pat) <- rownames(cyt)
WCT.v.nuc.con <- data.frame(c(WCT[1:4], nuc[1:4]))
rownames(WCT.v.nuc.con) <- rownames(cyt)

cyt.v.WCT.pat <- data.frame(c(cyt[5:10], WCT[5:10]))
rownames(cyt.v.WCT.pat) <- rownames(cyt)
cyt.v.WCT.con <- data.frame(c(cyt[1:4], WCT[1:4]))
rownames(cyt.v.WCT.con) <- rownames(cyt)

WCT_trunc.v.missense <- WCT[,5:10]
CYT_trunc.v.missense <- cyt[,5:10]
NUC_trunc.v.missense <- nuc[,5:10]

#Experimental input
exp <- WCT

#Rows must have at least 3 samples with scores of 10 or higher 
keep <- rowSums(exp>=10) >= 3
exp<-exp[keep,]

exp_info <- data.frame(condition = factor(c(rep("1", 4), rep("2", 6))))

# exp_info <- data.frame(condition = factor(c(rep("1", 6), rep("2", 6))),
                       # patientID = factor(c(1,2,3,4,5,6,1,2,3,4,5,6)))

#Create a coldata frame
coldata <- data.frame(row.names=colnames(exp), exp_info)
coldata 

#Create a DESeqDataSet object 
# dds <- DESeqDataSetFromMatrix(countData=exp, colData=coldata, design=~condition)
dds <- DESeqDataSetFromMatrix(countData=exp, colData=coldata, design=~condition)
dds

#Run DEseq2 pipleline
dds <- DESeq(dds)

res <- results(dds)
table(res$padj<0.05)
## Order by p-value
res <- res[order(res$pvalue), ]

# ## Merge with raw count data
# resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)

## Merge with normalised count data
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

#Convert ensembl IDs to HGNC symbols
ens_genes <- resdata_norm$Row.names
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=ens_genes,  mart=mart)


resdata_norm <- merge(resdata_norm, mart_back, by.x = "Row.names", by.y = "ensembl_gene_id")
genesort <- resdata_norm[order(resdata_norm$pvalue),]
genesort <- genesort[!duplicated(genesort$hgnc_symbol),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/PH_Fibroblasts/")
write.csv(genesort, "WCT_Alltranscripts_norm.csv", row.names = F)

Sig.padj <- subset(genesort, subset=(padj < 0.05))
Sig.padj.gene <- Sig.padj$hgnc_symbol
Sig.padj.gene <- Sig.padj.gene[!duplicated(Sig.padj.gene)]

Sig.padj.up <- subset(Sig.padj, subset=(log2FoldChange > 0))
Sig.padj.up.gene <- Sig.padj.up$hgnc_symbol
Sig.padj.up.gene <- Sig.padj.up.gene[!duplicated(Sig.padj.up.gene)]

Sig.padj.down <- subset(Sig.padj, subset=(log2FoldChange < 0))
Sig.padj.down.gene <- Sig.padj.down$hgnc_symbol
Sig.padj.down.gene <- Sig.padj.down.gene[!duplicated(Sig.padj.down.gene)]

write.table(Sig.padj.gene, "NUC_trunc.v.missense_sig_padj_genenames.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Sig.padj.up.gene, "NUC_trunc.v.missense_sig_padj_up_genenames.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Sig.padj.down.gene, "NUC_trunc.v.missense_sig_padj_down_genenames.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


##### TARDBP levels #####
setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/PH_Fibroblasts/")

cyt <- read.csv("Cytoplasm/Cytoplasm_norm.csv")
nuc <- read.csv("Nucleus/Nucleus_norm.csv")
wct <- read.csv("WCT/WCT_norm.csv")

rownames(cyt) <- cyt$hgnc_symbol
cytsample <- cyt[,8:17]

rownames(nuc) <- nuc$hgnc_symbol
nucsample <- nuc[,8:17]

rownames(wct) <- wct$hgnc_symbol
wctsample <- wct[,8:17]

cytTARDBP <- subset(cytsample, rownames(cytsample) == "TARDBP")
nucTARDBP <- subset(nucsample, rownames(nucsample) == "TARDBP")
WCTTARDBP <- subset(wctsample, rownames(wctsample) == "TARDBP")

TARDBP <- merge(cytTARDBP, nucTARDBP, by = 0)
TARDBP <- merge(TARDBP, WCTTARDBP, by.x = "Row.names", by.y = 0)
rownames(TARDBP) <- TARDBP$Row.names
TARDBP[,1] <- NULL

TARDBPv <- as.numeric(TARDBP[1,])
col = c(rep("chartreuse",4),
        rep("chartreuse4",6), 
        rep("dodgerblue",4),
        rep("dodgerblue4",6), 
        rep("gold2",4),
        rep("peru",6))

names <- colnames(TARDBP)
names <- sapply(strsplit(names,split= "_[aA-zZ]"),'[',1)

barplot(TARDBPv, 
        col = col,
        las = 2,
        cex.names = 1.1,
        names.arg = names, 
        ylab = "Normalised Counts")

legend("topright", 
       legend = c("Cytoplasm - Control", 
                  "Cytoplasm - Patient",
                  "Nucleus - Control",
                  "Nucleus - Patient", 
                  "WCT - Control",
                  "WCT - Patient"), 
       fill = c("chartreuse","chartreuse4","dodgerblue","dodgerblue4","gold2","peru"))



##### Results vs my DEGs #####

up_filter <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_up_genes.txt")
down_filter <-readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_down_genes.txt")
PPInetwork <-readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes.txt")


cytdeg <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/PH_Fibroblasts/Cytoplasm/Cyt_sig_padj_genenames.txt")
nucdeg <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/PH_Fibroblasts/Nucleus/Nuc_sig_padj_genenames.txt")
WCTdeg <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/PH_Fibroblasts/WCT/WCT_sig_padj_genenames.txt")

test <- WCTdeg

downoverlap <- Reduce(intersect, list(down_filter,test))
downoverlap
upoverlap <- Reduce(intersect, list(up_filter,test))
upoverlap
netoverlap <- Reduce(intersect, list(PPInetwork,test))
netoverlap


##### Results vs individual datasets #####

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian")
c9 <- read.csv("C9_unique.csv")
c91000 <- c9$Gene.Symbol[1:500]
sals <- read.csv("sals_unique.csv")
sals1000 <- sals$Gene.Symbol[1:500]
ftld <- read.csv("ftld_unique.csv")
ftld1000 <- ftld$Gene.Symbol[1:500]
VCP <- read.csv("VCP_unique.csv")
VCP1000 <- VCP$Gene.Symbol[1:500]

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_results_keepfiltering.csv")
pet1000 <- pet$hgnc_symbol[1:500]
rav <- read.csv("RAV_results_keepfiltering.csv")
rav1000 <- rav$hgnc_symbol[1:500]


C9overlap <- Reduce(intersect, list(c91000,test))
C9overlap

salsoverlap <- Reduce(intersect, list(sals1000,test))
salsoverlap

ftldoverlap <- Reduce(intersect, list(ftld1000,test))
ftldoverlap

vcpoverlap <- Reduce(intersect, list(vcp1000,test))
vcpoverlap

petoverlap <- Reduce(intersect, list(pet1000,test))
petoverlap

ravoverlap <- Reduce(intersect, list(rav1000,test))
ravoverlap


setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Collaborations/PH_Fibroblasts/Cyt.v.Nuc/")
down1 <- readLines("cyt_nuc_con_sig_padj_down_genenames.txt")
down2 <- readLines("nuc_cyt_pat_sig_padj_down_genenames.txt")

downoverlap <- Reduce(intersect, list(down1,down2))

up1 <- readLines("cyt_nuc_con_sig_padj_up_genenames.txt")
up2 <- readLines("nuc_cyt_pat_sig_padj_up_genenames.txt")

upoverlap <- Reduce(intersect, list(up1,up2))

resultsup <- subset(up2, !(up2 %in% upoverlap))
resultsdown <- subset(down2, !(down2 %in% downoverlap))

write.table(resultsup, "cyt.v.nuc_up_filtered.txt", row.names = F, quote = F)
write.table(resultsdown, "cyt.v.nuc_down_filtered.txt", row.names = F, quote = F)
