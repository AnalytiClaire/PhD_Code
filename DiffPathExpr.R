# ##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####

library (pathprint)
data(list = c("chipframe", "genesets","pathprint.Hs.gs","platform.thresholds", "pluripotents.frame"))
options(stringsAsFactors = FALSE)
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/NormalisedExpressionMatrices/") #set working directory to location of data

##### Run Single Chip Enrichment 2 #####
single.chip.enrichment2=function(exprs, geneset, statistic = "mean",normalizedScore = FALSE, progressBar = TRUE){
  # if (!(transformation %in% c("rank", "squared.rank", "log.rank"))) 
  #   stop("transformation should be rank, squared.rank or log.rank")
  # if (!(statistic %in% c("mean", "median"))) 
  #   stop("transformation should be mean or median")
  # if ((normalizedScore == TRUE & !(statistic == "mean"))) 
  #   stop("Parameteric normalization can only be used for statistic = mean")
  Ns <- ncol(exprs)
  Ng <- nrow(exprs)
  gene.names <- rownames(exprs)
  geneset.names <- names(geneset)
  # exprs <- apply(exprs, 2, rank, ties.method = "average")
  # if (transformation == "log.rank") {
  #   exprs <- log(exprs)
  # }
  # else if (transformation == "squared.rank") {
  #   exprs <- exprs^2
  # }
  if (progressBar == TRUE) {
    pb <- txtProgressBar(min = 0, max = length(geneset), 
                         style = 3)
  }
  score.matrix <- matrix(0, nrow = length(geneset), ncol = Ns)
  for (i in 1:length(geneset)) {
    overlap <- intersect(geneset[[i]], gene.names)
    if (length(overlap) == 0) {
      score.matrix[i, ] <- NA
    }
    else {
      if (statistic == "mean") {
        score.matrix[i, ] <- apply(exprs, 2, function(x) {
          mean(x[overlap])
        })
        # if (normalizedScore == TRUE) {  
        # 
        #   n <- length(overlap)
        #   if (transformation == "rank") {
        #     E.mean <- mean(1:Ng)
        #     E.sd <- ((sd(1:Ng)/(n^0.5))) * (((Ng - n)/(Ng - 
        #                                                  1))^0.5)
        #   }
        #   else if (transformation == "log.rank") {
        #     E.mean <- mean(log(1:Ng))
        #     E.sd <- ((sd(log(1:Ng))/(n^0.5))) * (((Ng - 
        #                                              n)/(Ng - 1))^0.5)
        #   }
        #   else if (transformation == "squared.rank") {
        #     E.mean <- mean((1:Ng)^2)
        #     E.sd <- ((sd((1:Ng)^2)/(n^0.5))) * (((Ng - 
        #                                             n)/(Ng - 1))^0.5)
        #   }
        #   score.matrix[i, ] <- sapply(score.matrix[i, 
        #                                            ], pnorm, mean = E.mean, sd = E.sd) - 0.5
        # }
      }
      else if (statistic == "median") {
        score.matrix[i, ] <- apply(exprs, 2, function(x) {
          median(x[overlap])
        })
      }
    }
    if (progressBar == TRUE) {
      setTxtProgressBar(pb, i)
    }
  }
  colnames(score.matrix) <- colnames(exprs)
  rownames(score.matrix) <- geneset.names
  return(score.matrix)
}




####C9_LCM ######
exp_C9 <- read.csv ("C9eset.csv", header=TRUE, row.names = 1)
C9ann <- customCDFAnn(exp_C9, chipframe$GPL570$ann)

EnrichmentC9 <- single.chip.enrichment2(exprs = C9ann,
                                         geneset = pathprint.Hs.gs,
                                         statistic = "mean",
                                         normalizedScore = FALSE,
                                         progressBar = TRUE
)

library(limma)

array_data.Control <- EnrichmentC9[,1:3]
array_data.Case <- EnrichmentC9[,4:11]

group = factor(c(colnames(array_data.Control),(colnames(array_data.Case)))) #substitute actual Case and Control datasets with array_data

Treat<-factor(rep(c("Control", "Patient"),c(3,8)), levels=c("Control", "Patient"))
design_df<-model.matrix(~Treat)
rownames(design_df)<-colnames(EnrichmentC9)
design_df



# contrasts_matrix = makeContrasts(TreatPatient,levels = design_df) #modify contrast levels as per group names in your dataset
# 
# fit<-lmFit(EnrichmentC9, design) #linear model fit
# fit2 = contrasts.fit(fit,contrasts_matrix)
# fit<-eBayes(fit) 
# result<-topTable(fit, coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 633) #"BH" adjust for multiple hypothesis testing


contrasts_matrix = makeContrasts(TreatPatient,levels = design_df) #modify contrast levels as per group names in your dataset
fit = lmFit(object = cbind(data.frame(array_data.Control,stringsAsFactors = F),data.frame(array_data.Case,stringsAsFactors = F)),design = design_df)
fit2 = contrasts.fit(fit,contrasts_matrix)
fit2=eBayes(fit2,trend = T)
diff_pathways=rownames(topTable(fit2,coef = 1,sort.by = "logFC",adjust.method = "BH",p.value = 0.1,number = 633))

result<-topTable(fit2, coef = 1,sort.by = "logFC",adjust.method = "BH", number = 633) #"BH" adjust for multiple hypothesis testing




# x <- subset(C9ann, rownames(C9ann) %in% pathprint.Hs.gs[[i]])
# C9ann <- C9ann[,4:11]

# ####sals_lcm###
exp_SALS <- read.csv ("sALSeset.csv", header=TRUE, row.names = 1)

# ####FTLD###
exp_FTLD <- read.csv ("FTLDeset.csv", header=TRUE, row.names = 1)

####VCP###
exp_VCP <- read.csv ("VCPeset.csv", header=TRUE, row.names = 1)

