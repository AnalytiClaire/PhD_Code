setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network")

Kegg_genes <- read.csv("KEGG_Genes.csv", na.strings = c("", "NA)"))
Kegg_genes<- as.list(Kegg_genes)
sub_list_all <- lapply(Kegg_genes, function(x) x[!is.na(x)])
sub_list <- sub_list_all[1:10]

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

write.csv(jac, "KEGG_jaccard_overlap_coefficient.csv")


#####
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian")

C9 <- read.csv("C9_unique.csv")
sALS <- read.csv("sALS_unique.csv")
FTLD <- read.csv("FTLD_unique.csv")
VCP <- read.csv("VCP_unique.csv")

pathway <- sub_list$Pathways.in.cancer

C9_merge <- subset(C9,C9$Gene.Symbol %in% pathway)







