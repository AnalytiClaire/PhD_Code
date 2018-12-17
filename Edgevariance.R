### Edge Variability of Nugget ###

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD")

edges <- read.csv("RMETHOD_samedir_mean.csv")
nugedge <- read.csv("NuggeteEdge.csv")
edges$swap <- paste(edges$Gene2,":",edges$Gene1, sep = "")

nugall <- subset(edges, edges$swap  %in% nugedge$Int)
nugdata <- nugall[,1:7]
rownames(nugdata) <- nugdata[,1]
nugdata[,1] <- NULL

plot.nug <- t(nugdata)

boxplot(plot.nug)


for (i in 1:ncol(plot.nug)){
  var.nug[i] <- var(plot.nug[,i])
}

nugall$variance <- var.nug
nugall <- nugall[order(nugall$variance),]

nugresult <- merge(nugedge, nugall, by.x = "Int", by.y = "swap")

write.csv(nugresult, "Edge_with_variance.csv", row.names = F)
