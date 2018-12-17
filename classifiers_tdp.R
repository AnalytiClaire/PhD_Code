# 
# #Install packages (if you don't already have them)
# install.packages(c('randomForest','ROCR'))
# install.packages(c('MASS','ggplot2','class'))
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("GEOquery")
# biocLite("Biobase")

# Loading up the packages
library(randomForest)
library(ROCR)
library(MASS)
library(ggplot2)
library(class)
library(Biobase)
library(GEOquery)
library(limma)

# Set working directory
setwd('/Users/clairegreen/Documents/PhD/ISB Summer Course/Day 3_Wednesday/Day 3_Afternoon')

##############################
### Load data for analysis ###
##############################
load("data/blood_miRNA_expression_LUSC_vs_Normal.rda")
View(train[,100:110])

SOD1 <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/non-TDP-43 Data Sets/JKSOD1_exprsgen.csv")

##############################
### Random Forest analysis ###
##############################
status <- vector
miRNA.rf = randomForest(as.factor(status) ~ .,data=subsetc9,importance=T,proximity=T)
#list variable importance score for the 109 miRNAs
miRNA.rf$importance[order(miRNA.rf$importance[,4]),]
# Plot variable importance
varImpPlot(miRNA.rf, type=2, n.var=10, main='Variable Importance for Top 10 Predictors\n(Mean decrease in Gini node impurity)')
# Get predictions for training data set
miRNA.predtrain <- predict(miRNA.rf,type="prob")
# Get predictions for test data set
miRNA.predtest <- predict(miRNA.rf,newdata=subsetsod1,type="prob")

######################################
### Calculate and plot ROC and AUC ###
######################################
pred.train = prediction(as.vector(miRNA.predtrain[,2]),as.vector(subsetc9[,c('status')]))
pred.train_auc = performance(pred.train, 'auc')
pred.train_rates = performance(pred.train, 'tpr','fpr')
plot(pred.train_rates, main='ROC miRNA Predictors', col='red', lwd=2)
pred.test = prediction(as.vector(miRNA.predtest[,2]),as.vector(subsetsod1[,c('status')]))
pred.test_auc = performance(pred.test, 'auc')
pred.test_rates = performance(pred.test, 'tpr','fpr')
plot(pred.test_rates, main='ROC miRNA Predictors', col='blue', lwd=2,add=T)
text(0.5,0.5,paste('AUC for train = ',format(pred.train_auc@y.values[[1]],digits=2,scientific=FALSE)),col="red")
text(0.5,0.4,paste('AUC for test = ',format(pred.test_auc@y.values[[1]],digits=2,scientific=FALSE)),col="blue")
grid()

##############################
### Select best n features ###
##############################
n <- 5
goodPredictors = rownames(miRNA.rf$importance)[order(miRNA.rf$importance[,4],decreasing=T)][1:n]

####################################
### Scatter plot of two features ###
####################################
f1 <- 2;
f2 <- 4;

matplot(train[train$status==0,goodPredictors[f1]],train[train$status==0,goodPredictors[f2]],
        main="Scatterplot of two best features",xlab=goodPredictors[f1],ylab=goodPredictors[f2],
        col='purple',pch=15)
matplot(train[train$status==1,goodPredictors[f1]],train[train$status==1,goodPredictors[f2]],
        main="Scatterplot of two best features",xlab=goodPredictors[f1],ylab=goodPredictors[f2],
        col='orange',pch=16,add=T)

