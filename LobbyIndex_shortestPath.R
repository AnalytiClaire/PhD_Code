setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/")

library(igraph)

Nodes <- readLines("TDP_fixed_NuggetGenes.txt")
Edges <- read.csv("CCM_Con_SOD1FUS_Removededge.csv")

g <- graph_from_data_frame(Edges, directed=FALSE, vertices=Nodes)
print(g, e=TRUE, v=TRUE)


#Lobby index
library(centiserve)
lob <- as.data.frame(lobby(g))

write.csv(lob, "TDPlobbyindex.csv")

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



setwd("/Users/clairegreen/Documents/PhD/Autoimmune/Results/TenRemove/")

library(igraph)

Nodes <- readLines("Knee_NugGenes.txt")
Edges <- read.csv("THISONE_AINug_Con_Tenpoint5_removededge.csv")

g <- graph_from_data_frame(Edges, directed=FALSE, vertices=Nodes)
print(g, e=TRUE, v=TRUE)


#Lobby index
library(centiserve)
lob <- as.data.frame(lobby(g))

write.csv(lob, "OAlobbyindex.csv")
