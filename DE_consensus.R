#####DIFFERENTIAL GENE EXPRESSION INTERSECT
#takes csv files of top X DE genes and identifies any consensus genes 

# setwd("/Users/clairegreen/Documents/PhD/TDP-43/non-TDP-43 Data Sets/JK2011_SOD1/")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/Thresholds/")
A <- read.csv(file =  "C9 _ap_6500.csv")

B <- read.csv(file =  "CH _ap_6500.csv")
 
C <- read.csv(file ="sALS _ap_6500.csv")

D <- read.csv(file ="FTLD _ap_6500.csv")

E <- read.csv(file = "VCP _ap_6500.csv")

F <- read.csv(file = "PETb_ap_6500.csv")

G <- read.csv(file = "RAV _ap_6500.csv")

# H <- read.csv(file =   "FSHD _ap_6500.csv")

A_DE <- A$Ensembl
B_DE <- B$Ensembl
C_DE <- C$Ensembl
D_DE <- D$Ensembl
E_DE <- E$Ensembl
F_DE <- F$Row.names
G_DE <- G$Row.names
# H_DE <- H$Gene.Symbol

overlap <- Reduce(intersect, list(B_DE, C_DE, D_DE, E_DE, F_DE, G_DE))
print(overlap)



setwd("/Users/clairegreen/Desktop/")
write.csv(overlap, file = "overlap_ens.csv")

