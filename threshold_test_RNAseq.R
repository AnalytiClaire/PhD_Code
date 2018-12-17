##RNA-Seq Gene Expression Analysis using Limma##

library(limma)
library(edgeR)
library(biomaRt)
library(tictoc)

analysis.name<-"Pet" #Label analysis
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/Ravits/")
# Counts <- read.table(file = 'GSE67196_Petrucelli2015_ALS_genes.rawcount.txt', header = TRUE)
# 
# write.csv(x = Counts, file = "counts_petrucelli.csv")

Counts <- read.csv(file = "ravitscombined.counts.csv", header = TRUE)
# Counts <- as.matrix(Counts[-28])
m <- as.matrix(Counts[,-1])
rownames(m) <- Counts[,1]
Counts <- m


# num_overlap <- matrix(data=NA)

#####
##RUN FROM HERE

deglen <- list()
genlen <- list()
tic()
for (i in 1:500){
  
  tic()
  Countnum <- Counts[rowSums(Counts) >i,]
  #DGElist
  dge <- DGEList(counts=Countnum)
  dge <- calcNormFactors(dge)
  #Design
  Treat<-factor(rep(c("Control", "Patient"),c(8,13)), levels=c("Control", "Patient"))
  design<-model.matrix(~Treat)
  rownames(design)<-colnames(Countnum)
  
  #Voom transformation
  v <- voom(dge,design,plot=FALSE)
  
  #Limma fitting
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(Countnum)) #"BH" adjust for multiple hypothesis testing

  deglen[[i]] <- length(which(result$P.Value < 0.05))
  genlen[[i]] <- length(result$logFC)
  print(paste("Loop finished for iteration: ", i))
  toc()
}
toc()
# genlen_df <- t(as.data.frame(genlen))
# deglen_df <- t(as.data.frame(deglen))



plot(genlen_df, deglen_df)

genlen <- as.numeric(genlen)
deglen <- as.numeric(deglen)
comb <- cbind(genlen, deglen)
comb <- transform.default(comb)
# comb <- transform(comb, deglen / genlen )
comb[, "ratio"] <- comb[, "deglen"] / comb[, "genlen"]

plot(genlen, comb$ratio)
