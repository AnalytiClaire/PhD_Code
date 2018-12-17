setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/")
net <- read.csv("Merged Network_DEG_C9_VCP_GRN default node.csv")

C9 <- subset(net, net$C9 == "true")
VCP <- subset(net, net$VCP == "true")
GRN <- subset(net, net$GRN == "true")

c9name <- as.character(C9$name)
write.table(c9name, "C9_genenames.txt", row.names = F, col.names = F, quote = F)

vcpname <- as.character(VCP$name)
write.table(vcpname, "VCP_genenames.txt", row.names = F, col.names = F, quote = F)

grnname <- as.character(GRN$name)
write.table(grnname, "GRN_genenames.txt", row.names = F, col.names = F, quote = F)

Reduce(intersect, list(vcpname, grnname))

C9only <- data.frame(subset(c9name, !(c9name %in% c(vcpname, grnname))))
C9only$C9only <- "true"
names(C9only) <-c("Genename", "C9only")

VCPonly <- data.frame(subset(vcpname, !(vcpname %in% c(c9name, grnname))))
VCPonly$VCPonly <- "true"
names(VCPonly) <-c("Genename", "VCPonly")

GRNonly <- data.frame(subset(grnname, !(grnname %in% c(vcpname, c9name))))
GRNonly$GRNonly <- "true"
names(GRNonly) <-c("Genename", "GRNonly")

net <- merge(net, C9only, by.x = "shared.name", by.y = "Genename")


write.csv(C9only, "C9only_genenames.csv",row.names = F, col.names = F, quote = F)
write.csv(VCPonly, "VCPonly_genenames.csv",row.names = F, col.names = F, quote = F)
write.csv(GRNonly, "GRNonly_genenames.csv",row.names = F, col.names = F, quote = F)
