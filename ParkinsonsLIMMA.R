##Differential Expression of Genes##

library(affy)
library(tcltk)
library(widgetTools)
library(DynDoc)
library(tools)
library(Biobase)
library(tkWidgets)
library(plyr)
###########################
### Sporadic Parkinson's ##
###########################

##### LEW GSE19587 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE19587_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"LEW" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(4,6)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)

##### MOR.FC GSE8397 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE8397_Cortex/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MOR.FC" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(3,5)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)



##### MID3 GSE20168 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20168_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MID3" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(15,14)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)





##### MID4 GSE20291 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20291_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MID4" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(20,15)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)





##### MID1 GSE20141 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20141_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"DAV" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(8,10)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)


##### FFR GSE7621 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE7621_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"FFR" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(9,16)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)
##### MID2 GSE20292 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20292_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MID" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(18,11)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)


##### DIJ GSE49036 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE49036_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"DIJ" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(8,15)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)


##### MOR.SN GSE8397 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE8397_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MOR" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
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
Treat <- factor(rep(c("Control", "Patient"),c(15,24)), levels=c("Control", "Patient"))
Tissue <- factor(c(rep(c("SN", "LN"),c(7,9)),rep(c("SN", "LN"),c(8,15))), levels=c("SN", "LN"))

design<-model.matrix(~Treat+Tissue)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)



##### MOUSE ######
##### CON GSE4758 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE4758/")
annotation <- read.csv("Mouse430_2.na36.ortholog.csv")
annotation_HG <- subset(annotation, annotation$Ortholog.Array == "HG-U133_Plus_2" )

annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
annotation_U133aplus2<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)


annotation_HG$Ortholog.Probe.Set <- tolower(annotation_HG$Ortholog.Probe.Set)
annotation_merge <- merge(annotation_HG, annotation_U133aplus2, by.x = "Ortholog.Probe.Set", by.y = "Probe.Set.ID")
annotation_merge[,8] <- NULL
annotation_merge[,14] <- NULL
# write.table(annotation_merge, "Mouse430_2.na36.ortholog_merge_HG_U133A_Plus_2.txt", row.names = F)

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"CON" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation <- annotation_merge
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
Treat <- factor(rep(c("Control", "Patient"),c(6,6)), levels=c("Control", "Patient"))
# Tissue <- factor(c(rep(c("SN", "LN"),c(7,9)),rep(c("SN", "LN"),c(8,15))), levels=c("SN", "LN"))

design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- as.character(result$Gene.Symbol)
result$Ensembl <- as.character(result$Ensembl)
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,8]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)




###########################
#### LRRK2 Parkinson's ####
###########################
#### BOT GSE34516 #########
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/")
Data <- read.csv("GSE34516.RMA-GENE-EXTENDED-EC-hg19-na36_LRRK2only.csv", row.names = 1)

genename <- as.data.frame(Data$GeneSymbol)
rownames(genename) <- rownames(Data)
colnames(genename) <- "Gene.Symbol"
#Data is already normalised using RMA-sketch
analysis.name<-"BOT" #Label analysis
expressionMatrix<-Data[,1:6] #takes expression from normalised expression set


Treat<-factor(rep(c("Control", "Patient"),c(4,2)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))


result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by = 0) #merge values into one array
result<-merge(result, genename, by.x = "Row.names", by.y = 0)
result<-subset(result, subset=(Gene.Symbol !="---")) #if no gene symbol, discount


setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,16]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

#### BOT2 GSE23290 #########
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/")
Data <- read.csv("GSE23290.RMA-GENE-EXTENDED-EC-hg19-na36_LRRK2only.csv", row.names = 1)

genename <- as.data.frame(Data$GeneSymbol)
rownames(genename) <- rownames(Data)
colnames(genename) <- "Gene.Symbol"
#Data is already normalised using RMA-sketch
analysis.name<-"BOT2" #Label analysis
expressionMatrix<-Data[,1:8] #takes expression from normalised expression set


Treat<-factor(rep(c("Control", "Patient"),c(5,3)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))


result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by = 0) #merge values into one array
result<-merge(result, genename, by.x = "Row.names", by.y = 0)
result<-subset(result, subset=(Gene.Symbol !="---")) #if no gene symbol, discount


setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,18]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)






###########################
#### Parkinson's Blood ####
###########################
##### AMA GSE99039 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE99039_keepfiles/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"AMA" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(183,169)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

uniqueresult<-subset(uniqueresult, Gene.Symbol!="") #removes any probes for which there are no gene symbols
uniqueresult<-subset(uniqueresult, subset=(countP>20)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)

##### RON GSE72267 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE72267_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"RON" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(19,40)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)









##### AMA ATP13A2 #########
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE99039_keepfiles/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"AMA_ATP13A2" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(183,5)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

uniqueresult<-subset(uniqueresult, Gene.Symbol!="") #removes any probes for which there are no gene symbols
uniqueresult<-subset(uniqueresult, subset=(countP>20)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)

##### AMA PINK1 ##########
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE99039_keepfiles/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"AMA_PINK1" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(183,12)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

uniqueresult<-subset(uniqueresult, Gene.Symbol!="") #removes any probes for which there are no gene symbols
uniqueresult<-subset(uniqueresult, subset=(countP>20)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)

##### AMA PARKIN ##########
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE99039_keepfiles/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"AMA_PARKIN" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


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

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
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
Treat<-factor(rep(c("Control", "Patient"),c(183,13)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

uniqueresult<-subset(uniqueresult, Gene.Symbol!="") #removes any probes for which there are no gene symbols
uniqueresult<-subset(uniqueresult, subset=(countP>20)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)
###########################
# Parkinson's Fibroblasts #
###########################
######### LRRK2Fib ########
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/")
Data <- read.csv("LRRK2_PARKIN_Fibroblasts-RMA-GENE-EXTENDED-EC-hg19-na36.csv", row.names = 1)

genename <- as.data.frame(Data$Gene.Symbol)
rownames(genename) <- rownames(Data)
colnames(genename) <- "Gene.Symbol"
#Data is already normalised using RMA-sketch
analysis.name<-"LRRK2Fib" #Label analysis
expressionMatrix<-Data[,1:9] #takes expression from normalised expression set


Treat<-factor(rep(c("Control", "Patient"),c(6,3)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))


result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by = 0) #merge values into one array
result<-merge(result, genename, by.x = "Row.names", by.y = 0)
result<-subset(result, subset=(Gene.Symbol !="---")) #if no gene symbol, discount


setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,16]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

######### PARK2Fib ########
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/")
Data <- read.csv("LRRK2_PARKIN_Fibroblasts-RMA-GENE-EXTENDED-EC-hg19-na36.csv", row.names = 1)

genename <- as.data.frame(Data$GeneSymbol)
rownames(genename) <- rownames(Data)
colnames(genename) <- "Gene.Symbol"
#Data is already normalised using RMA-sketch
analysis.name<-"PARK2Fib" #Label analysis
expressionMatrix<-Data[,c(1:3,10:12)] #takes expression from normalised expression set


Treat<-factor(rep(c("Control", "Patient"),c(3,3)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))


result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by = 0) #merge values into one array
result<-merge(result, genename, by.x = "Row.names", by.y = 0)
result<-subset(result, subset=(Gene.Symbol !="---")) #if no gene symbol, discount


setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,16]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

