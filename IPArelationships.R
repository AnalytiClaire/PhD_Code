#IPA relationships
library(dplyr)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD")

IPA <- read.csv("IPA_Relationships.csv", stringsAsFactors = F)
nug <- readLines("NuggetGenes.txt")

sub.either <- subset(IPA, IPA$From.Molecule.s. %in% nug | IPA$To.Molecule.s. %in% nug)

sub.both <- subset(IPA, IPA$From.Molecule.s. %in% nug & IPA$To.Molecule.s. %in% nug)

rmrows <- which(sub.both[,1] == sub.both[,3])
rmrows <- paste(rmrows, collapse=", ")

sub.both.unique <- sub.both[-c(rmrows),]

write.csv(sub.both, "IPA_interactions_nuggetonly.csv")
write.csv(sub.both.unique, "IPA_interactions_Unique.csv")
