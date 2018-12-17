setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/")
library(igraph)

Nodes <- readLines("PD_NuggetGenes_conoverlapremoved.txt")
Edges <- read.csv("PDNugget_conremove_edge.csv")

g <- graph_from_data_frame(Edges, directed=FALSE, vertices=Nodes)
print(g, e=TRUE, v=TRUE)


#Lobby index
library(centiserve)
lob <- as.data.frame(lobby(g))

write.csv(lob, "PDlobbyindex.csv")

#

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network")
Nodes <- read.csv("DEG_PPI_default_node.csv")
Edges <- read.csv("DEG_PPI_Network default edge.csv")
Edges <- Edges[,c(6,7,5)]

g <- graph_from_data_frame(d = Edges, directed = F, vertices = Nodes)

test <- list()

for (i in 1:length(Nodes)){
  <- g[[i]]
}


gdist <- distances(g)

SP <- shortest_paths(graph = g, from = 13, to = 2072)

ASP <- all_shortest_paths(graph = g, from = 13, to = 2072)
ASP <- as.data.frame(ASP)
