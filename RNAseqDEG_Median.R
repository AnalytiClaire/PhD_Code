##RNA-Seq Gene Expression Analysis using Limma##

analysis.name<-"RAV" #Label analysis
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/Ravits/")
# Counts <- read.table(file = 'GSE67196_Petrucelli2015_ALS_genes.rawcount.txt', header = TRUE)
# 
# write.csv(x = Counts, file = "counts_petrucelli.csv")

Counts <- read.csv(file = "ravitsannotated_combined.counts.csv", header = TRUE)

Counts[Counts == 0] <- NA
# Counts[Counts<30] <- NA
Counts <- na.omit(Counts)
rownames(Counts)<-Counts[,1]
Counts[,1] <- NULL

# Counts<-subset(Counts, subset=(GeneID !="NA")) #if no gene symbol, discount

# Countszero <-subset(Counts, subset=(row !=0))
# Countszero <- apply(Counts, 1, function(row) all(row !="NA"))
# Counts <- Counts[Countszero,]

library(limma)
library(edgeR)

#remove gene symbols from end column
Countnum <- Counts[,1:26]
# Counts <- data.matrix(Counts)

# Countnum <- read.csv(file = "pet.counts.clean.csv")

#DGElist
dge <- DGEList(counts=Countnum)
dge <- calcNormFactors(dge)
dge

#Design
Treat<-factor(rep(c("Control", "Patient"),c(8,13)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(Countnum)
design

#Voom transformation
v <- voom(dge,design,plot=FALSE)

#Limma fitting
fit <- lmFit(v,design)
fit <- eBayes(fit)
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(Countnum)) #"BH" adjust for multiple hypothesis testing
result <- merge(result, Counts, by="row.names", all=TRUE)
#result <- result[,1:7]

#Count tables from bcbio have ensembl gene IDs. This must be annotated with HGNC symbols

#Download the HGNC symbols and gene IDs using a vector containing the IDs from results
library(biomaRt)
genes <- as.vector(result[,1])
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)

# library(org.Hs.eg.db)
# library(GeneNetworkBuilder)

#Merge the tables using ensembl ID
result <- merge(result, mart_back, by.x = "Row.names", by.y = "ensembl_gene_id")
# result[,1] <- NULL




#### Take median value for gene duplicates ###########
result3 <- ddply(result,"hgnc_symbol", numcolwise(median, (result$adj.P.Val)))
#result3 <- aggregate(result, by=list("Gene.Symbol"), FUN=median)

genesort <- result[order(result$P.Value),]

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/MedianGenes")
write.csv(genesort, file=paste(analysis.name, "ENSrankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)
```

microarray, filename "Wenbin_DE_Gene.R"

```{r, eval=FALSE}
##Differential Expression of Genes##

library(affy)
library(tcltk)
library(widgetTools)
library(DynDoc)
library(tools)
library(Biobase)
library(tkWidgets)
library(plyr)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/FTLD/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
#celfiles<-basename(celfiles)
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"FTLD" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set

#write.csv(dataMatrixAll, file = "eset.GRN.csv")

#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation file

###USING BIOMART
# library (biomaRt)
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org") 
# x <- rownames(dataMatrixAll) #create vector containing probe IDs
# mart_attribute <- listAttributes(mart)
# annotation <- getBM(attributes=c("affy_hg_u133a_2", "hgnc_symbol", "description"), 
#                    filters = "affy_hg_u133a_2", values = x, mart = mart)
# annotation<-subset(annotation, subset=(hgnc_symbol !="")) #if no gene symbol, discount

# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
#annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35.annot.txt"
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot.txt"
#annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(8,16)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all


result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")

# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
#write.csv(result, file=paste(analysis.name, "result.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe

#### Take median value for gene duplicates ###########
result2 <- ddply(result,"Gene.Symbol", numcolwise(median, (result$adj.P.Val)))
#result3 <- aggregate(result, by=list("Gene.Symbol"), FUN=median)

genesort <- result2[order(result2$adj.P.Val),]

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/MedianGenes")
write.csv(genesort, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)