library(tictoc)
library(gdata)

load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint_biomart.RData")
setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")
pathways <- read.csv("KEGG_Pathways_uniquemutation.csv", header = F)
pathways <- pathways$V1

# sub_list <- subset(pathprint_list, names(pathprint_list) %in% pathways)
sub_list <- read.csv(file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/GRN/GRN_KEGG.csv", na.strings = c("", "NA)"))
sub_list <- as.list(sub_list)
sub_list <- lapply(sub_list, function(x) x[!is.na(x)])

#Make function for finding Jaccard coefficient
jaccard<-function(A,B){
  jc=length(intersect(A,B))/length(union(A,B))
  return(jc)
}

#Create empty matrix 
jc_mat <- matrix(0, length(sub_list), length(sub_list))

tic()
for (i in 1:length(sub_list)){
    jc_mat[i,] <-unname(unlist(lapply(sub_list,function(x)jaccard(A = sub_list[[i]],B = x))))
  }

rownames(jc_mat) <- colnames(jc_mat) <- names(sub_list)
toc()


##Only take upper triangle without diagonal (all comparisons are currently doubled)
ptri <- jc_mat
ptri[lower.tri(ptri, diag = TRUE)] <- NA

# Turn into vector
vec <- unmatrix(ptri)
# Remove NA values
vec <- na.omit(vec)
jac <- as.data.frame(vec)

write.csv(jac, "GRN_jaccard_overlap_coefficient.csv")
