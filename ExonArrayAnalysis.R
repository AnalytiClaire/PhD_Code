library(oligo)
library(affy)
library(tcltk)
library(widgetTools)
library(DynDoc)
library(tools)
library(Biobase)
library(tkWidgets)
library(plyr)

#Method taken from:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3220870/ 

##### Command Line #####
#Download Human Exon 1.0 ST Array Analysis (zip, 131 MB) and Human Exon 1.0 ST Array Probeset, and Meta Probeset Files, core, full, extended and comprehensive hg18 (zip, 13 MB)
#from affymetrix.com

#Download APT from https://www.thermofisher.com/uk/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html

#Exon based
#/Users/clairegreen/Documents/PhD/Parkinsons/apt-1.21.0-x86_64-apple-yosemite/bin/apt-probeset-summarize -a rma-sketch -p /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.pgf -c /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.clf -s /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.core.ps -qc-probesets /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.qcc -o OUT_EXON *.CEL

#Gene Based
/Users/clairegreen/Documents/PhD/Parkinsons/apt-1.21.0-x86_64-apple-yosemite/bin/apt-probeset-summarize -a rma-sketch -p /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.pgf -c /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.clf -m /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2-dt1-hg18-ps/HuEx-1_0-st-v2.r2.dt1.hg18.core.mps -qc-probesets /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.qcc -o OUT_GENE *.CEL

#####  #####

setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE34516_RAW/OUT_GENE")
qc <- read.table("rma-sketch.report.txt", sep="\t", header=T)
plot(1:10, qc$pm_mean, ylim=c(0,1000), xlab="Array", ylab="Signal Intensity", main="Average Raw Intensity Signal")
plot(1:10, qc$all_probeset_mad_residual_mean, ylim=c(0,0.3), xlab="Array", ylab="Mean absolute deviation", main="Deviation of Residuals from Median")


d.gene <- read.table("rma-sketch.summary.txt", sep="\t", header=T, row.names=1)
d.t <- dist(t(d.gene))
plot(hclust(d.t), main="Hierarchical clustering", labels=c(rep("Control", 4), rep("LRRK2", 2), rep("Sporadic", 4)))
plot(density(d.gene[,1]), main="Distribution of RMA-normalised intensities", xlab="RMA normalised intensity")
for(i in 2:ncol(d.gene)) {lines (density(d.gene[,i]))}

##### Compute P Values

/Users/clairegreen/Documents/PhD/Parkinsons/apt-1.21.0-x86_64-apple-yosemite/bin/apt-probeset-summarize -a dabg -p /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.pgf -c /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.clf -b /Users/clairegreen/Documents/PhD/Parkinsons/HuEx-1_0-st-v2-r2/HuEx-1_0-st-v2.r2.antigenomic.bgp -o ./OUT_DABG *.CEL


##
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE34516_RAW/OUT_DABG/")

dabg <- read.table("dabg.summary.txt", sep="\t", header=T, row.names=1)
dim(dabg) # 1411399 10
dabg.core <- subset(dabg, rownames(dabg) %in% rownames(d.gene))
dim(dabg.core) # 287329 10




setwd("/Users/clairegreen/Documents/PhD/LRRK2_PARKIN_Project/CEL files/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-read.celfiles(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"LEW" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set