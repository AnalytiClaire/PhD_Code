### RNA-seq using DEseq2 on Count Matrix

library(DESeq2)

#Read in count matrix
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h")

Counts <- read.csv(file = "combined.counts.csv", header = TRUE)

rownames(Counts)<-Counts[,1]
Counts[,1] <- NULL

## Divide into different analyses
TRL <- Counts[,1:6]
WCT <- Counts[,13:18]
CYT <- Counts[,19:24]
GFP_L <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","GRASPS.TRL.GFPLow.1",
                   "GRASPS.TRL.GFPLow.2","GRASPS.TRL.GFPLow.3")]
GFP_H <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","GRASPS.TRL.GFPHigh.1",
                   "GRASPS.TRL.GFPHigh.2","GRASPS.TRL.GFPHigh.3")]
GRASPS_VS_WCT <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","WCT.Control.1",
                           "WCT.Control.2","WCT.Control.3")]
GRASPS_VS_CYT <- Counts[,c("GRASPS.TRL.Control.1","GRASPS.TRL.Control.2","GRASPS.TRL.Control.3","CytTrc.Control.1",
                           "CytTrc.Control.2","CytTrc.Control.3")]
WCT_VS_CYT <- Counts[,c("WCT.Control.1","WCT.Control.2","WCT.Control.3","CytTrc.Control.1",
                        "CytTrc.Control.2","CytTrc.Control.3")]

#Select experiment
analysis.name <- "TRL"
countdata <- TRL

countdata[rowSums(countdata) == 0,] <- NA
countdata <- na.omit(countdata)

#Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

#Assign condition
condition <- factor(c(rep("Condition1", 3), rep("Condition2", 3)))

#Create a coldata frame
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata
#Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

#Run DEseq2 pipleline
dds <- DESeq(dds)

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/DEseq2/")
# # Plot dispersions
# png("qc-dispersions.png", 1000, 1000, pointsize=20)
# plotDispEsts(dds, main="Dispersion plot")
# dev.off()
# 
# # Regularized log transformation for clustering/heatmaps, etc
# rld <- rlogTransformation(dds)
# head(assay(rld))
# hist(assay(rld))
# 
# # Colors for plots below
# ## Ugly:
# ## (mycols <- 1:length(unique(condition)))
# ## Use RColorBrewer, better
# library(RColorBrewer)
# (mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
# 
# # Sample distance heatmap
# sampleDists <- as.matrix(dist(t(assay(rld))))
# library(gplots)
# png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
# heatmap.2(as.matrix(sampleDists), key=F, trace="none",
#           col=colorpanel(100, "black", "white"),
#           ColSideColors=mycols[condition], RowSideColors=mycols[condition],
#           margin=c(15, 15), main="Sample Distance Matrix")
# dev.off()
# 
# # Principal components analysis
# ## Could do with built-in DESeq2 function:
# ## DESeq2::plotPCA(rld, intgroup="condition")
# ## I like mine better:
# rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
#   require(genefilter)
#   require(calibrate)
#   require(RColorBrewer)
#   rv = rowVars(assay(rld))
#   select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
#   pca = prcomp(t(assay(rld)[select, ]))
#   fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
#   if (is.null(colors)) {
#     if (nlevels(fac) >= 3) {
#       colors = brewer.pal(nlevels(fac), "Paired")
#     }   else {
#       colors = c("black", "red")
#     }
#   }
#   pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
#   pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
#   pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
#   pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
#   plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
#   with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
#   legend(legendpos, legend=levels(fac), col=colors, pch=20)
#   #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
#   #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
#   #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
# }
# png("qc-pca.png", 1000, 1000, pointsize=20)
# rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
# dev.off()
# 


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata_norm <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

genes <- as.vector(resdata$Gene)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)
result <- merge(resdata, mart_back, by.x = "Gene", by.y = "ensembl_gene_id")
result_med <- ddply(result,"hgnc_symbol", numcolwise(median, (result$padj)))
#result3 <- aggregate(result, by=list("Gene.Symbol"), FUN=median)

genesort <- result_med[order(result_med$padj),]

## Write results
write.csv(genesort, file=paste(analysis.name, "_diffexpr-results.csv", sep=""), row.names=FALSE, quote = FALSE)


CG_TRL_deseq <- genesort

### 

CG_TRL_sig_deseq <- subset(CG_TRL_deseq, subset=(padj < 0.05))
CG_WCT_gene_deseq <- CG_TRL_sig_deseq$hgnc_symbol

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Guillaume_HEK48h/DEseq2/")
write.table(CG_WCT_gene_deseq, "CG_WCT_gene_deseq.txt", row.names=FALSE, quote = FALSE)



overlap <- Reduce(intersect, list(CG_WCT_gene, common_DEGs))
print(overlap)



write.table(CG_TRL_gene, "CG_TRL_Gene_DESeq2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)










## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()


