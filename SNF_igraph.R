### SNF igraph ###

load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/SimilarityNetworkFusion/Microarray_RNAseq/WGCNA_DS0_Min8/DS0_Min8.RData")

test <- cyt
out <- cyt
col <- geneInfo

for (i in 1:nrow(test)){
  x <- test[i,]
  y1 <- merge(x,col, by.x = "gene1", by.y = "Gene")
  y2 <- merge(x,col, by.x = "gene2", by.y = "Gene")
  out[i,4] <- y1$module
  out[i,5] <- y2$module
  if (grepl("^.+(00)$", i)){
    message(paste(i), " rows completed")
  } else if (i == nrow(cyt)) {
    message("All rows completed")
  }
}

colors <- unique(dynamicColors)
den_table <- list()
PC_table <- list()
AC_table <- list()

for (i in 1:length(colors)){
  PC <- length(which(col == colors[i]))
  PC <- (PC*(PC - 1))/2
  mod <- subset(out, out$V4 == colors[i] & out$V5 == colors[i])
  AC <- nrow(mod)
  den <- AC/PC
  
  den_table[i] <- den
  PC_table[i] <- PC
  AC_table[i] <- AC
  
  names(den_table)[i] <- colors[i]
  names(PC_table)[i] <- colors[i]
  names(AC_table)[i] <- colors[i]
}



#Whole network density (density = number of actual connections/number of potential connections)
PC.net <- nrow(col) #how many nodes in the network
PC.net <- (PC.net*(PC.net - 1))/2 #how many potential connections
AC.net <- nrow(cyt) #how many actual connections
den.net <- AC.net/PC.net #calculate density
den.net


#Extra-cluster density
PC.ex <- nrow(col) #number of total nodes
PC.ex <- (PC.ex*(PC.ex - 1))/2 - Reduce("+", den_table) #number of possible connections - possible intra-cluster connections
AC.ex <- nrow(out) - nrow(out[out$V4 == out$V5,])
den.ex <- AC.ex/PC.ex 
den.ex

#Intra-cluster density
PC.in <- Reduce("+", PC_table)
AC.in <- Reduce("+", AC_table)
den.in <- AC.in/PC.in 
den.in

#Average cluster density
Reduce("+", den_table)/length(colors)


