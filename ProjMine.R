
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/ProjMine/NewGenes_Oct2018/")

dataset2015 <- read.csv("Patients-2015.vcf.avinput.hg19_multianno.txtclaire.muts.als.csv")
dataset2016 <- read.csv("Patients-2016.vcf.avinput.hg19_multianno.txtclaire.muts.als.csv")
dataset2017 <- read.csv("Patients-2017.vcf.avinput.hg19_multianno.txtclaire.muts.als.csv")
datasetmisc <- read.csv("Patients.vcf.avinput.hg19_multianno.txtclaire.muts.als.csv")

data_all <- rbind(dataset2015, dataset2016, dataset2017, datasetmisc)

#Merge and then replace far right columns with new column information
# data <- merge(dataset2015,dataset2016, by = "Start", all = T)
# write.csv(data, "test.csv")
# 
# data2 <- read.csv("test.csv")
# data2 <- merge(data2,dataset2017, by = "Start", all = T)
# write.csv(data2, "test2.csv")
# 
# data3 <- read.csv("test2.csv")
# data3 <- merge(data3,datasetmisc, by = "Start", all = T)
# write.csv(data3, "test3.csv")

# data <- read.csv("All_ProjMine.csv")


##### Dataset 2015 #####
data_filter <- dataset2015
# data_filter <- data[which(data$DANN_score > 0.96),]
data_filter[data_filter == "./.:.:.:.:.:."] <- "NA"
data_filter[data_filter == "./.:.:.:.:."] <- "NA"

data_search <- data_filter
rownames(data_search) <- data_search$Start
data_search <- data_search[,10:ncol(data)]
data_search <- as.matrix(data_search)

mutcount <- apply(data_search, 1, function(x) length(which(!is.na(x))))
mutcount <- as.data.frame(mutcount, row.names = names(mutcount))

Pat2015 <- data_filter[1:9]
Pat2015 <- merge(Pat2015, mutcount, by.x = "Start", by.y = 0)
write.csv(Pat2015, "Pat2015.csv")

##### Dataset 2016  #####
data <- dataset2016
data_filter <- data[which(data$DANN_score > 0.96),]
data_filter[data_filter == "./.:.:.:.:.:.:."] <- "NA"
data_filter[data_filter == "./.:.:.:.:.:."] <- "NA"

data_search <- data_filter
rownames(data_search) <- data_search$Start
data_search <- data_search[,10:ncol(data)]
data_search <- as.matrix(data_search)

mutcount <- apply(data_search, 1, function(x) length(which(!is.na(x))))
mutcount <- as.data.frame(mutcount, row.names = names(mutcount))

Pat2016 <- data_filter[1:9]
Pat2016 <- merge(Pat2016, mutcount, by.x = "Start", by.y = 0)
write.csv(Pat2016, "Pat2016.csv")

##### Dataset 2017  #####
data <- dataset2017
data_filter <- data[which(data$DANN_score > 0.96),]
data_filter[data_filter == "./.:.:.:.:.:."] <- "NA"
data_filter[data_filter == "./.:.:.:.:."] <- "NA"

data_search <- data_filter
rownames(data_search) <- data_search$Start
data_search <- data_search[,10:ncol(data)]
data_search <- as.matrix(data_search)

mutcount <- apply(data_search, 1, function(x) length(which(!is.na(x))))
mutcount <- as.data.frame(mutcount, row.names = names(mutcount))

Pat2017 <- data_filter[1:9]
Pat2017 <- merge(Pat2017, mutcount, by.x = "Start", by.y = 0)
write.csv(Pat2017, "Pat2017.csv")

##### Dataset Misc  #####
data <- datasetmisc
data_filter <- data[which(data$DANN_score > 0.96),]
data_filter[data_filter == "./.:.:.:.:.:."] <- "NA"
data_filter[data_filter == "./.:.:.:.:."] <- "NA"

data_search <- data_filter
rownames(data_search) <- data_search$Start
data_search <- data_search[,10:ncol(data)]
data_search <- as.matrix(data_search)

mutcount <- apply(data_search, 1, function(x) length(which(!is.na(x))))
mutcount <- as.data.frame(mutcount, row.names = names(mutcount))

PatMisc <- data_filter[1:9]
PatMisc <- merge(PatMisc, mutcount, by.x = "Start", by.y = 0)
write.csv(PatMisc, "PatMisc.csv")



controlmisc <- read.csv("Cntrl.vcf.avinput.hg19_multianno.txtclaire.muts.als.csv")
control2016 <- read.csv("Cntrl-2016.vcf.avinput.hg19_multianno.txtclaire.muts.als.csv")
control2017 <- read.csv("Cntrl-2017.vcf.avinput.hg19_multianno.txtclaire.muts.als.csv")


##### Control Misc  #####
data <- controlmisc
data_filter <- data[which(data$DANN_score > 0.96),]
data_filter[data_filter == "./.:.:.:.:.:."] <- "NA"
data_filter[data_filter == "./.:.:.:.:."] <- "NA"

data_search <- data_filter
rownames(data_search) <- data_search$Start
data_search <- data_search[,10:ncol(data)]
data_search <- as.matrix(data_search)

mutcount <- apply(data_search, 1, function(x) length(which(!is.na(x))))
mutcount <- as.data.frame(mutcount, row.names = names(mutcount))

conmisc <- data_filter[1:9]
conmisc <- merge(conmisc, mutcount, by.x = "Start", by.y = 0)
write.csv(conmisc, "controlmisc.csv")

##### Control 2016  #####
data <- control2016
data_filter <- data[which(data$DANN_score > 0.96),]
data_filter[data_filter == "./.:.:.:.:.:.:."] <- "NA"

data_search <- data_filter
rownames(data_search) <- data_search$Start
data_search <- data_search[,10:ncol(data)]
data_search <- as.matrix(data_search)

mutcount <- apply(data_search, 1, function(x) length(which(!is.na(x))))
mutcount <- as.data.frame(mutcount, row.names = names(mutcount))

con2016 <- data_filter[1:9]
con2016 <- merge(con2016, mutcount, by.x = "Start", by.y = 0)
write.csv(con2016, "control2016.csv")

##### Control 2017  #####
data <- control2017
data_filter <- data[which(data$DANN_score > 0.96),]
data_filter[data_filter == "./.:.:.:.:.:."] <- "NA"

data_search <- data_filter
rownames(data_search) <- data_search$Start
data_search <- data_search[,10:ncol(data)]
data_search <- as.matrix(data_search)

mutcount <- apply(data_search, 1, function(x) length(which(!is.na(x))))
mutcount <- as.data.frame(mutcount, row.names = names(mutcount))

con2017 <- data_filter[1:9]
con2017 <- merge(con2017, mutcount, by.x = "Start", by.y = 0)
write.csv(con2017, "control2017.csv")


#Patient Merge

#Open each test file and rank gene name in alphabetical order. For NAs, cut and paste the "y" info into the "x" columns. 
#Remove the extra Y information KEEPING THE MUTATION COUNTS COLUMNS

result <- merge(Pat2015,Pat2016, by = "Start", all = T)
write.csv(result, "test.csv", row.names = F)

data2 <- read.csv("test.csv")
data2 <- merge(data2,Pat2017, by = "Start", all = T)
write.csv(data2, "test2.csv", row.names = F)

data3 <- read.csv("test2.csv", stringsAsFactors = F)
data3 <- merge(data3,PatMisc, by = "Start", all = T)
write.csv(data3, "test3.csv", row.names = F)

#Patient Merge
result <- merge(con2016,con2017, by = "Start", all = T)
write.csv(result, "test.csv", row.names = F)

data2 <- read.csv("test.csv")
data2 <- merge(data2,conmisc, by = "Start", all = T)
write.csv(data2, "test2.csv", row.names = F)



#Merge with topological information

patresult <- read.csv("PM_Pat_all.csv")
conresult <- read.csv("PM_Con_All.csv")

nodeinfo <- read.csv("Nugget_default_node.csv")

Patmerge <- merge(patresult,nodeinfo, by.x = "Gene.refGene.x", by.y = "name", all = T)
Conmerge <- merge(conresult, nodeinfo,by.x = "Gene.refGene.x", by.y = "name", all = T )

write.csv(Patmerge, "Patient+nodeinfo.csv", row.names = F)
write.csv(Conmerge, "Controls+nodeinfo.csv", row.names = F)


#Remove any common rows

Pat_Start <- Patmerge$Start
Con_Start <- Conmerge$Start

CommonStart <- intersect(Pat_Start, Con_Start)

PatUnique <- Patmerge[!(Patmerge$Start %in% CommonStart),]

PatUnique$PatFreq <- sapply(PatUnique$Total, function(x) x/1152)

PatUnique$ExAC_ALL.x <- replace(PatUnique$ExAC_ALL.x,is.na(PatUnique$ExAC_ALL.x),0)

PatFinal <- PatUnique[!(PatUnique$ExAC_ALL.x > PatUnique$PatFreq),]
# PatFinal <- PatFinal[!(is.na(PatFinal$ExAC_ALL.x)),]
PatFinal$ExAC_ALL.x <- replace(PatFinal$ExAC_ALL.x, ,values = "NA")


write.csv(PatFinal, "PatFinal.csv", row.names = F)

##### Merge with edge number

final <- read.csv("PatFinal.csv")

#aggregate number of patients with ANY SNP in gene
final_ag <- aggregate(final$Total, by=list(Category=final$Gene.refGene.x), FUN=sum)

#aggregate number of SNPs per gene
final$one <- 1
final_SNPsperGene <- aggregate(final$one, by=list(Category=final$Gene.refGene.x), FUN=sum)


#read in edges
edges <- read.csv("Nugget_default_node.csv")

#Merge with edge info (patients)
final_ag_edge <- merge(final_ag, edges, by.x = "Category", by.y = "name")
rownames(final_ag_edge) <- final_ag_edge$Category
final_ag_edge$Category <- NULL

#Merge with edge info (SNPS)
final_SNP_edge <- merge(final_SNPsperGene, edges, by.x = "Category", by.y = "name")
rownames(final_SNP_edge) <- final_SNP_edge$Category
final_SNP_edge$Category <- NULL

#Test correlation
cor.test(final_ag_edge$x, final_ag_edge$NeighborhoodConnectivity)
cor.test(final_SNP_edge$x, final_SNP_edge$BetweennessCentrality)

######

ID <- read.csv("PatientId.csv", stringsAsFactors = F)
sev <- read.csv("Claire.network.survival.csv")


rownumbers <- which(ID$PatID2015 %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatID2015[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatID2015[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatID2016 %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatID2016[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatID2016[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatID2016b %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatID2016b[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatID2016b[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatID2017 %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatID2017[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatID2017[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatID2017b %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatID2017b[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatID2017b[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatID2017c %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatID2017c[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatID2017c[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatIDmisc %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatIDmisc[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatIDmisc[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatIDmiscb %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatIDmiscb[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatIDmiscb[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatIDmiscc %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatIDmiscc[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatIDmiscc[rownumbers[i]] <- x
}

rownumbers <- which(ID$PatIDmiscd %in% sev$PatientID)
for (i in 1:length(rownumbers)){
  y <- ID$PatIDmiscd[rownumbers[i]]
  sevrow <- which(sev$PatientID == y)
  x <- as.numeric(sev[sevrow,2])
  ID$PatIDmiscd[rownumbers[i]] <- x
}

colnames(ID) <- c("Start", "Ref", "Alt", "mutcount2015", "mutcount2016", "mutcount2017", "mutcountmisc", "Sev2015", "Sev2016",
                  "Sev2016b", "Sev2017", "Sev2017b", "Sev2017c", "SevMisc", "SevMiscb", "SevMiscc", "SevMiscd")
write.csv(ID, "IDtoSeverity.csv", row.names = F)

ID <-read.csv("IDtoSeverity.csv")
nodeinfo <- read.csv("Patient+nodeinfo.csv")

Allinfo <- merge(nodeinfo, ID, by = "Start", all = T)
Allinfo[,28:33] <- NULL

Severity <- as.matrix(Allinfo[,28:37])
Severity<- apply(Severity, 2, as.numeric)

Allinfo$MeanSeverity <- x
x <- rowMeans(Severity, na.rm=TRUE)

cor.test(Allinfo_noNA$MeanSeverity, Allinfo$BetweennessCentrality)


Allinfo_noNA <- Allinfo[!is.na(Allinfo$MeanSeverity),]
MeanperGene_SNP <- aggregate(Allinfo_noNA$MeanSeverity, by=list(Category=Allinfo_noNA$Gene.refGene.x), FUN=median)

MeanSev_edge <- merge(MeanperGene, edges, by.x = "Category", by.y = "name")

cor.test(MeanSev_edge$BetweennessCentrality, MeanSev_edge$x)

####
nona <- read.csv("AllInfo_noNA.csv")

final_unique <- subset(final, !(final$Start %in% nona$Start))

write.csv(final_unique, "claireisstupid.csv", row.names = F)
