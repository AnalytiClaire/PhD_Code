library(psych)
library(tictoc)
library(gdata)

#Load data
uniqueresult <- read.csv("C9result.csv")
rownames(uniqueresult) <- uniqueresult[,1]
uniqueresult[,1] <- NULL


##transpose
CorExprMat <- t(uniqueresult)


#Generate correlation and p values
tic()
cortest <- corr.test(CorExprMat, use = "pairwise", method = "pearson", adjust = "fdr") #must use transposed matrix (genes are colnames)
toc()

#Extract R values
cortestoutput <- cortest$r
corRadj <- cortestoutput
corRadj[lower.tri(corRadj, diag = TRUE)] <- NA

#Turn into vector
corRadj <- as.matrix(corRadj)
corRvec <- unmatrix(corRadj)
#Remove NA values
corRvec <- na.omit(corRvec)
corRvec <- as.data.frame(corRvec)


#Extract P values
cortestpadjust <- cortest$p
corPadj <- cortestpadjust
corPadj[lower.tri(corPadj, diag = TRUE)] <- NA

#Turn into vector
corPadj <- as.matrix(corPadj)
corPvec <- unmatrix(corPadj)
#Remove NA values
corPvec <- na.omit(corPvec)
corPvec <- as.data.frame(corPvec)

#merge into one data frame
CorData <- merge(corRvec, corPvec, by.x = "row.names", by.y = "row.names")

#Select significant results
sigoutput <- subset(CorData, CorData$corPvec < 0.05)
write.csv(sigoutput, file = "sigcor.csv")
