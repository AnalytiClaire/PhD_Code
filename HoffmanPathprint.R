##### Using pathprint to identify common pathways across multiple TDP-43 pathology-containing data sets ####

library (pathprint)
options(stringsAsFactors = FALSE)


setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HoffmanMuscle/") #set working directory to location of data
exprs.matrix <- read.csv ("eset.hoffman.csv", header=TRUE) #assign the .csv file to a variable, column headers are true
row.names (exprs.matrix) <- exprs.matrix[,1] #specify that first column contains gene names
exprs.matrix<- exprs.matrix[,2:(ncol(exprs.matrix))] #specify that all other columns are gene expression data

pathprint.result <- exprs2fingerprint(exprs.matrix, platform = "GPL96", species="human", progressBar=T)
vec.data <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Diff.pathways <- diffPathways(pathprint.result, vec.data, 0.5)

print(Diff.pathways)
