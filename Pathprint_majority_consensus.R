load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprinted.RData")
thresh = (100/50)
#Select patient columns
pat_C9 <- pathprint_C9[,4:11]
pat_sals <- pathprint_sals[,4:10]
pat_ftld <- pathprint_ftld[,9:24]
pat_vcp <- pathprint_vcp[,4:10]
#Combine
pat_all_pp <- cbind(pat_C9, pat_sals,pat_ftld, pat_vcp)

keeppat <- rowSums(pat_all_pp==1) >= ncol(pat_all_pp)/thresh
pat_up<-pat_all_pp[keeppat,]
keeppat <- rowSums(pat_all_pp==0) >= ncol(pat_all_pp)/thresh
pat_same<-pat_all_pp[keeppat,]
keeppat <- rowSums(pat_all_pp==-1) >= ncol(pat_all_pp)/thresh
pat_down<-pat_all_pp[keeppat,]


#Select control columns
con_C9 <- pathprint_C9[,1:3]
con_sals <- pathprint_sals[,1:3]
con_ftld <- pathprint_ftld[,1:8]
con_vcp <- pathprint_vcp[,1:3]
#Combine
con_all_pp <- cbind(con_C9, con_sals,con_ftld, con_vcp)

keepcon <- rowSums(con_all_pp==1) >= ncol(con_all_pp)/thresh
con_up<-con_all_pp[keepcon,]
keepcon <- rowSums(con_all_pp==0) >= ncol(con_all_pp)/thresh
con_same<-con_all_pp[keepcon,]
keepcon <- rowSums(con_all_pp==-1) >= ncol(con_all_pp)/thresh
con_down<-con_all_pp[keepcon,]


#Create summary column for pathways based on their grouping
pat1 <- data.frame(row.names = rownames(pat_up))
pat1$Values <- 1
pat0 <- data.frame(row.names = rownames(pat_same))
pat0$Values <- 0
patneg1 <- data.frame(row.names = rownames(pat_down))
patneg1$Values <- -1
#Merge dataframes together
pat_all <- rbind(pat1, pat0, patneg1)

#Create summary column for pathways based on their grouping
con1 <- data.frame(row.names = rownames(con_up))
con1$Values <- 1
con0 <- data.frame(row.names = rownames(con_same))
con0$Values <- 0
conneg1 <- data.frame(row.names = rownames(con_down))
conneg1$Values <- -1
#Merge dataframes together
con_all <- rbind(con1, con0, conneg1)


pp_all <- merge(pat_all, con_all, by = 0)


diff_consensus_MAJ <- subset(pp_all, !(pp_all$Values.x == pp_all$Values.y))

down_path <- subset(diff_consensus_MAJ, diff_consensus_MAJ$Values.x < diff_consensus_MAJ$Values.y)
down_pathways <- down_path$Row.names
up_path <- subset(diff_consensus_MAJ, diff_consensus_MAJ$Values.x > diff_consensus_MAJ$Values.y)
up_pathways <- up_path$Row.names

dyspathways_MAJ <- c(down_pathways, up_pathways)
