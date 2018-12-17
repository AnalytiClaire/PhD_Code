setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15")
library(DCGL)
exprs <- read.csv("C9rankeduniqueresult.csv")

## divide exprs into two parts corresponding to condition 1
##(exprs.1) and condition 2 (exprs.2) respectively
rownames(exprs)<-exprs$Gene.Symbol
exprs<-exprs[,49:59]
exprs.1<-exprs[1:3]
exprs.2<-exprs[4:11]
DCp.res<-DCp(exprs.1,exprs.2,
             r.method='pearson',
             link.method='qth',cutoff=0.25,N=0)
iDCe.res<-DCe(exprs.1,exprs.2,
             link.method='qth',
             cutoff=0.25,
             nbins=10,p=0.1)
## combine two Differential Co-expression Analysis results
DCsum.res<-DCsum(DCp.res,DCe.res,
                 DCpcutoff=0.25,DCecutoff=0.4)
DCsum.res$DCGs[1:3,]
DCsum.res$DCLs[1:3,]
## sort out differentially regulated genes and differentially regulated links
data(tf2target) ## TF-to-target relationships
DRsort.res<-DRsort(DCsum.res$DCGs,DCsum.res$DCLs,tf2target,expGenes)
## or
DRsort.res<-DRsort(DCe.res$DCGs,DCe.res$DCLs,tf2target,expGenes)
## plot differentially regulated links
DRplot.res<-DRplot(DCsum.res$DCGs,DCsum.res$DCLs,
                   tf2target,
                   expGenes,
                   4 ASC
                   type='TF_bridged_DCL',
                   vsize=5,asize=0.25,lcex=0.3,ewidth=1,
                   figname=c('TF2target_DCL.pdf','TF_bridged_DCL.pdf'))
## rank regulators by TED or TDD
DRrank.res<-DRrank(DCsum.res$DCGs,DCsum.res$DCLs,
                   tf2target,
                   expGenes,
                   rank.method=c('TED','TDD')[1],
                   Nperm=0)
## rank regulators by RIF\
data(exprs_design)
RIF.res<-RIF(exprs,exprs.1,exprs.2,
             tf2target,
             exprs_design,
             p.value=0.05)
ASC Identi