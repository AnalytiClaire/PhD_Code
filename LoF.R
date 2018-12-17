setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/Prioritisation/")

nug <- readLines("../TDP_fixed_NuggetGenes.txt")
lof1 <- read.csv("All_data_ALS_1.csv")
lof2 <- read.csv("All_data_ALS_2.csv")


lof <- rbind(lof1, lof2)
lof_rare <- subset(lof, lof$is_lof == 1 & lof$aaf_exac_all < 0.01 & lof$aaf_1kg_all < 0.01)
lof_nug <- subset(lof_rare, lof_rare$gene %in% nug)

write.csv(lof_nug, "Nug_LOF_ProjectMinE.csv", row.names = F)


lof_novel <- subset(lof, lof$is_lof == 1 & lof$in_1kg == 0 & lof$in_exac == 0)
lof_nug <- subset(lof_novel, lof_novel$gene %in% nug)
write.csv(lof_nug, "Nug_novelLOF_ProjectMinE.csv", row.names = F)
