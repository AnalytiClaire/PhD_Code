library(pathprint)

load("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/allpathprint.RData")
names <- read.table("~/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Consensus_fingerprint/dysregulatedpathways.txt")
names <- as.vector(names$V1)

#Select patient columns
pat_C9 <- pathprint_C9[,4:11]
pat_C9 <- subset(pat_C9, row.names(pat_C9) %in% names)
con_C9 <- pathprint_C9[,1:3]
con_C9 <- subset(con_C9, row.names(con_C9) %in% names)

pat_sals <- pathprint_sals[,4:10]
pat_sals <- subset(pat_sals, row.names(pat_sals) %in% names)
con_sals <- pathprint_sals[,1:3]
con_sals <- subset(con_sals, row.names(con_sals) %in% names)

pat_ftld <- pathprint_ftld[,9:24]
pat_ftld <- subset(pat_ftld, row.names(pat_ftld) %in% names)
con_ftld <- pathprint_ftld[,1:8]
con_ftld <- subset(con_ftld, row.names(con_ftld) %in% names)

pat_vcp <- pathprint_vcp[,4:10]
pat_vcp <- subset(pat_vcp, row.names(pat_vcp) %in% names)
con_vcp <- pathprint_vcp[,1:3]
con_vcp <- subset(con_vcp, row.names(con_vcp) %in% names)

pat_FUS <- pathprint_fus[,4:6]
pat_FUS <- subset(pat_FUS, row.names(pat_FUS) %in% names)
con_FUS <- pathprint_fus[,1:3]
con_FUS <- subset(con_FUS, row.names(con_FUS) %in% names)

pat_SOD1 <- pathprint_sod1[,8:10]
pat_SOD1 <- subset(pat_SOD1, row.names(pat_SOD1) %in% names)
con_SOD1 <- pathprint_sod1[,1:7]
con_SOD1 <- subset(con_SOD1, row.names(con_SOD1) %in% names)

pat_CBFTLD <- pathprint_cbftld[,8:16]
pat_CBFTLD <- subset(pat_CBFTLD, row.names(pat_CBFTLD) %in% names)
con_CBFTLD <- pathprint_cbftld[,1:7]
con_CBFTLD <- subset(con_CBFTLD, row.names(con_CBFTLD) %in% names)



#Combine
all_pp <- cbind(pat_C9, pat_sals,pat_ftld, pat_vcp,pat_FUS, pat_SOD1, pat_CBFTLD,
                    con_C9, con_sals,con_ftld, con_vcp,con_FUS, con_SOD1, con_CBFTLD)



pat_TDP_pp <- cbind(pat_C9, pat_sals,pat_ftld, pat_vcp)
con_TDP_pp <- cbind(con_C9, con_sals, con_ftld, con_vcp)
pat_nonTDP_pp <- cbind(pat_SOD1, pat_FUS)
con_nonTDP_pp <- cbind(con_SOD1, con_FUS)

pat_TDPCB_pp<- cbind(pat_C9, pat_sals,pat_ftld, pat_vcp, pat_CBFTLD)
pat_nonTDPCB_pp<- cbind( pat_FUS, pat_SOD1,pat_CBFTLD)

colnames(pat_all_pp) <- c("C9_1", "C9_2", "C9_3",
                          "C9_4", "C9_5", "C9_6",
                          "C9_7", "C9_8",
                          "sALS_1", "sALS_2", "sALS_3", "sALS_4",
                          "sALS_5", "sALS_6", "sALS_7",
                          "FTLD_FC_1","FTLD_FC_2","FTLD_FC_3","FTLD_FC_4",
                          "FTLD_FC_5","FTLD_FC_6","FTLD_FC_7","FTLD_FC_8",
                          "FTLD_FC_9","FTLD_FC_10","FTLD_FC_11","FTLD_FC_12",
                          "FTLD_FC_13","FTLD_FC_14","FTLD_FC_15","FTLD_FC_16",
                          "VCP_1","VCP_2","VCP_3","VCP_4","VCP_5","VCP_6",
                          "VCP_7",
                          "FUS_1", "FUS_2", "FUS_3",
                          "SOD1_1", "SOD1_2","SOD1_3",
                          "FTLD_CB_1","FTLD_CB_2","FTLD_CB_3","FTLD_CB_4",
                          "FTLD_CB_5","FTLD_CB_6","FTLD_CB_7","FTLD_CB_8",
                          "FTLD_CB_9")

colnames(pat_some_pp) <- c("C9_1", "C9_2", "C9_3",
                          "C9_4", "C9_5", "C9_6",
                          "C9_7", "C9_8",
                          "sALS_1", "sALS_2", "sALS_3", "sALS_4",
                          "sALS_5", "sALS_6", "sALS_7",
                          "FTLD_FC_1","FTLD_FC_2","FTLD_FC_3","FTLD_FC_4",
                          "FTLD_FC_5","FTLD_FC_6","FTLD_FC_7","FTLD_FC_8",
                          "FTLD_FC_9","FTLD_FC_10","FTLD_FC_11","FTLD_FC_12",
                          "FTLD_FC_13","FTLD_FC_14","FTLD_FC_15","FTLD_FC_16",
                          "VCP_1","VCP_2","VCP_3","VCP_4","VCP_5","VCP_6",
                          "VCP_7",
                          "FUS_1", "FUS_2", "FUS_3",
                          "SOD1_1", "SOD1_2","SOD1_3")

colnames(pat_TDP_pp) <- c("C9_1", "C9_2", "C9_3",
                           "C9_4", "C9_5", "C9_6",
                           "C9_7", "C9_8",
                           "sALS_1", "sALS_2", "sALS_3", "sALS_4",
                           "sALS_5", "sALS_6", "sALS_7",
                           "FTLD_FC_1","FTLD_FC_2","FTLD_FC_3","FTLD_FC_4",
                           "FTLD_FC_5","FTLD_FC_6","FTLD_FC_7","FTLD_FC_8",
                           "FTLD_FC_9","FTLD_FC_10","FTLD_FC_11","FTLD_FC_12",
                           "FTLD_FC_13","FTLD_FC_14","FTLD_FC_15","FTLD_FC_16",
                           "VCP_1","VCP_2","VCP_3","VCP_4","VCP_5","VCP_6",
                           "VCP_7")
colnames(pat_TDPCB_pp) <- c("C9_1", "C9_2", "C9_3",
                          "C9_4", "C9_5", "C9_6",
                          "C9_7", "C9_8",
                          "sALS_1", "sALS_2", "sALS_3", "sALS_4",
                          "sALS_5", "sALS_6", "sALS_7",
                          "FTLD_FC_1","FTLD_FC_2","FTLD_FC_3","FTLD_FC_4",
                          "FTLD_FC_5","FTLD_FC_6","FTLD_FC_7","FTLD_FC_8",
                          "FTLD_FC_9","FTLD_FC_10","FTLD_FC_11","FTLD_FC_12",
                          "FTLD_FC_13","FTLD_FC_14","FTLD_FC_15","FTLD_FC_16",
                          "VCP_1","VCP_2","VCP_3","VCP_4","VCP_5","VCP_6",
                          "VCP_7",
                          "FTLD_CB_1","FTLD_CB_2","FTLD_CB_3","FTLD_CB_4",
                          "FTLD_CB_5","FTLD_CB_6","FTLD_CB_7","FTLD_CB_8",
                          "FTLD_CB_9")

colnames(pat_nonTDP_pp) <- c("FUS_1", "FUS_2", "FUS_3",
                           "SOD1_1", "SOD1_2","SOD1_3")

colnames(pat_nonTDPCB_pp) <- c("FUS_1", "FUS_2", "FUS_3",
                             "SOD1_1", "SOD1_2","SOD1_3",
                             "FTLD_CB_1","FTLD_CB_2","FTLD_CB_3","FTLD_CB_4",
                             "FTLD_CB_5","FTLD_CB_6","FTLD_CB_7","FTLD_CB_8",
                             "FTLD_CB_9")




thresh=0.5
#Run consensus
pat_consensus <- data.frame(consensusFingerprint(pat_TDP_pp, thresh))
con_consensus <- data.frame(consensusFingerprint(con_TDP_pp, thresh))
conDis <- consensusDistance(con_consensus, con_CBFTLD)
conDis$samples <-rownames(conDis)
# hist(conDis$distance,conDis$samples)
dotchart(conDis$distance, labels = row.names(conDis), 
         cex = 0.5, main = "FTLD cerebellum compared to non-TDP consensus")
head(conDis)


col = c("cyan", "white", "red")

Condition = HeatmapAnnotation(df = data.frame(type1 = c(rep("Patient", 56), rep("Control", 40))), 
                              col = list(type1 = c("Patient" =  "red", "Control" = "blue")), show_legend = TRUE)
Disease = HeatmapAnnotation(df = data.frame(Condition = c(rep("C9_Pat", 8), 
                                                          rep("sALS_Pat", 7), 
                                                          rep("FTLD_Pat", 16),
                                                          rep("VCP_Pat", 7),
                                                          rep("FUS_Pat", 3),
                                                          rep("SOD1_Pat", 3),
                                                          rep("CBFTLD_Pat",9),
                                                          rep("C9_Con", 3),
                                                          rep("sALS_Con", 3),
                                                          rep("FTLD_Con", 8),
                                                          rep("VCP_Con", 3),
                                                          rep("FUS_Con", 3),
                                                          rep("SOD1_Con", 7),
                                                          rep("CBFTLD_Con",7))), 
                            col = list(Condition = c("C9_Pat" =  "olivedrab1", 
                                                     "C9_Con" = "olivedrab3", 
                                                     "FTLD_Pat" = "skyblue1", 
                                                     "FTLD_Con" = "skyblue3",
                                                     "sALS_Pat" = "gold", 
                                                     "sALS_Con" = "gold3",
                                                     "VCP_Pat" = "orchid2", 
                                                     "VCP_Con" = "orchid3",
                                                     "FUS_Pat" = "azure3", 
                                                     "FUS_Con" = "azure4",
                                                     "SOD1_Pat" = "firebrick1", 
                                                     "SOD1_Con" = "firebrick3",
                                                     "CBFTLD_Pat" = "royalblue1",
                                                     "CBFTLD_Con"= "navy")), show_legend = TRUE)

Heatmap(all_pp,
        name = "Pathprint Score",
        col = col,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        clustering_method_rows = "complete",
        clustering_distance_rows = "euclidean",
        show_row_names = TRUE,
        show_column_names = FALSE,
        show_row_dend = TRUE,
        top_annotation = Condition, 
        column_names_max_height = unit(60, "mm"))


###################################
##### Patient FTLD Cerebellum #####
###################################

#Vs Con Cerebellum
consensus <- data.frame(consensusFingerprint(con_CBFTLD, thresh))
patCB_conCB <- consensusDistance(consensus, pat_CBFTLD)

#Vs TDP+pat
consensus <- data.frame(consensusFingerprint(pat_TDP_pp, thresh))
patCB_patTDP <- consensusDistance(consensus, pat_CBFTLD)

#Vs TDP-pat
consensus <- data.frame(consensusFingerprint(pat_nonTDP_pp, thresh))
patCB_patnonTDP <- consensusDistance(consensus, pat_CBFTLD)

#Vs TDP+con
consensus <- data.frame(consensusFingerprint(con_TDP_pp, thresh))
patCB_conTDP <- consensusDistance(consensus, pat_CBFTLD)

#Vs TDP-con
consensus <- data.frame(consensusFingerprint(con_nonTDP_pp, thresh))
patCB_connonTDP <- consensusDistance(consensus, pat_CBFTLD)

                                     
###################################
##### Control FTLD Cerebellum #####
###################################

#Vs TDP+pat
consensus <- data.frame(consensusFingerprint(pat_TDP_pp, thresh))
conCB_patTDP <- consensusDistance(consensus, con_CBFTLD)

#Vs TDP-pat
consensus <- data.frame(consensusFingerprint(pat_nonTDP_pp, thresh))
conCB_patnonTDP <- consensusDistance(consensus, con_CBFTLD)

#Vs TDP+con
consensus <- data.frame(consensusFingerprint(con_TDP_pp, thresh))
conCB_conTDP <- consensusDistance(consensus, con_CBFTLD)

#Vs TDP-con
consensus <- data.frame(consensusFingerprint(con_nonTDP_pp, thresh))
conCB_connonTDP <- consensusDistance(consensus, con_CBFTLD)

###### EXTRACT DISTANCES #######
patCB_conCB_dis <-patCB_conCB$distance
patCB_patTDP_dis <-patCB_patTDP$distance
patCB_patnonTDP_dis <-patCB_patnonTDP$distance
patCB_conTDP_dis <-patCB_conTDP$distance
patCB_connonTDP_dis <-patCB_connonTDP$distance
conCB_patTDP_dis <-conCB_patTDP$distance
conCB_patnonTDP_dis <-conCB_patnonTDP$distance
conCB_conTDP_dis <-conCB_conTDP$distance
conCB_connonTDP_dis <-conCB_connonTDP$distance

#### Take mean distance ####
patCB_conCB_m <-mean(patCB_conCB_dis)
patCB_patTDP_m <-mean(patCB_patTDP_dis)
patCB_patnonTDP_m <-mean(patCB_patnonTDP_dis)
patCB_conTDP_m <-mean(patCB_conTDP_dis)
patCB_connonTDP_m <-mean(patCB_connonTDP_dis)
conCB_patTDP_m <-mean(conCB_patTDP_dis)
conCB_patnonTDP_m <-mean(conCB_patnonTDP_dis)
conCB_conTDP_m <-mean(conCB_conTDP_dis)
conCB_connonTDP_m <-mean(conCB_connonTDP_dis)

patCB_conCB_m
patCB_patTDP_m
patCB_patnonTDP_m
patCB_conTDP_m
patCB_connonTDP_m
conCB_patTDP_m
conCB_patnonTDP_m
conCB_conTDP_m
conCB_connonTDP_m

  ###### Extract PValues #######
patCB_conCB_p <-patCB_conCB$pvalue
patCB_patTDP_p <-patCB_patTDP$pvalu
patCB_patnonTDP_p <-patCB_patnonTDP$pvalue
patCB_conTDP_p <-patCB_conTDP$pvalue
patCB_connonTDP_p <-patCB_connonTDP$pvalue
conCB_patTDP_p <-conCB_patTDP$pvalue
conCB_patnonTDP_p <-conCB_patnonTDP$pvalue
conCB_conTDP_p <-conCB_conTDP$pvalue
conCB_connonTDP_p <-conCB_connonTDP$pvalue

#### Take mean distance ####
patCB_conCB_mp <-mean(patCB_conCB_p)
patCB_patTDP_mp <-mean(patCB_patTDP_p)
patCB_patnonTDP_mp <-mean(patCB_patnonTDP_p)
patCB_conTDP_mp <-mean(patCB_conTDP_p)
patCB_connonTDP_mp <-mean(patCB_connonTDP_p)
conCB_patTDP_mp <-mean(conCB_patTDP_p)
conCB_patnonTDP_mp <-mean(conCB_patnonTDP_p)
conCB_conTDP_mp <-mean(conCB_conTDP_p)
conCB_connonTDP_mp <-mean(conCB_connonTDP_p)

patCB_conCB_mp
patCB_patTDP_mp
patCB_patnonTDP_mp
patCB_conTDP_mp
patCB_connonTDP_mp
conCB_patTDP_mp
conCB_patnonTDP_mp
conCB_conTDP_mp
conCB_connonTDP_mp
