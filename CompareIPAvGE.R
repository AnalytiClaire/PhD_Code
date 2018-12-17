setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/GeneXplain")

pw <- read.csv("Pathways4PCxN.csv")
pw$TermIPA <- sub("$", " (KEGG)", pw$TermIPA) #This means "Add at the end ($) of the string the following

pw$TermIPA <- toupper(pw$TermIPA)
pw$TermGE <- toupper(pw$TermGE)

write.table(intersect(pw$PCxN, pw$TermIPA), "PCxN_IPA.txt", quote = F, row.names = F, col.names = F)
write.table(intersect(pw$PCxN, pw$TermGE),"PCxN_GE.txt", quote = F, row.names = F, col.names = F)
write.table(intersect(pw$TermGE, pw$TermIPA),"GE_IPA.txt", quote = F, row.names = F, col.names = F)

