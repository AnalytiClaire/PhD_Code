### PATIENT VARIANCE ####
options(scipen=999)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
rownames(C9) <- C9$Gene.Symbol
C9 <- C9[,12:19]

sals <- read.csv("sals_unique.csv")
rownames(sals) <- sals$Gene.Symbol
sals <- sals[,12:18]

ftld <- read.csv("ftld_unique.csv")
rownames(ftld) <- ftld$Gene.Symbol
ftld <- ftld[,17:32]

vcp <- read.csv("vcp_unique.csv")
rownames(vcp) <- vcp$Gene.Symbol
vcp <- vcp[,12:18]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_results_keepfiltering.csv")
rownames(pet) <- pet$hgnc_symbol
pet <- pet[,19:35]

rav <- read.csv("RAV_results_keepfiltering.csv")
rownames(rav) <- rav$hgnc_symbol
rav <- rav[,18:30]

## Calculate rowise variance
c9var <- apply(C9, 1, var)
salsvar <- apply(sals, 1, var)
ftldvar <- apply(ftld, 1, var)
vcpvar <- apply(vcp, 1, var)
petvar <- apply(pet, 1, var)
ravvar <- apply(rav, 1, var)

varC9 <- data.frame(row.names = rownames(C9), 
                    variance <- c9var)
varsals <- data.frame(row.names = rownames(sals), 
                    variance <- salsvar)
varftld <- data.frame(row.names = rownames(ftld), 
                    variance <- ftldvar)
varvcp <- data.frame(row.names = rownames(vcp), 
                    variance <- vcpvar)
varpet <- data.frame(row.names = rownames(pet), 
                     variance <- petvar)
varrav <- data.frame(row.names = rownames(rav), 
                     variance <- ravvar)

names <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/intersect_up_1.txt")
names <- names$V1 

up_varC9 <- subset(varC9, rownames(varC9) %in% names)
up_varsals<- subset(varsals, rownames(varsals) %in% names)
up_varftld <- subset(varftld, rownames(varftld) %in% names)
up_varvcp <- subset(varvcp, rownames(varvcp) %in% names)
up_varpet <- subset(varpet, rownames(varpet) %in% names)
up_varrav <- subset(varrav, rownames(varrav) %in% names)

up_varC9$gene_symbol <- rownames(up_varC9)
up_varC9 <- up_varC9[order(up_varC9$variance....c9var),]
up_varC9$rank <- 1:nrow(up_varC9)

up_varsals$gene_symbol <- rownames(up_varsals)
up_varsals <- up_varsals[order(up_varsals$variance....salsvar),]
up_varsals$rank <- 1:nrow(up_varsals)

up_varftld$gene_symbol <- rownames(up_varftld)
up_varftld <- up_varftld[order(up_varftld$variance....ftldvar),]
up_varftld$rank <- 1:nrow(up_varftld)

up_varvcp$gene_symbol <- rownames(up_varvcp)
up_varvcp <- up_varvcp[order(up_varvcp$variance....vcpvar),]
up_varvcp$rank <- 1:nrow(up_varvcp)

up_varpet$gene_symbol <- rownames(up_varpet)
up_varpet <- up_varpet[order(up_varpet$variance....petvar),]
up_varpet$rank <- 1:nrow(up_varpet)

up_varrav$gene_symbol <- rownames(up_varrav)
up_varrav <- up_varrav[order(up_varrav$variance....ravvar),]
up_varrav$rank <- 1:nrow(up_varrav)

C9rank <- up_varC9[order(up_varC9$gene_symbol),]
salsrank <- up_varsals[order(up_varsals$gene_symbol),]
ftldrank <- up_varftld[order(up_varftld$gene_symbol),]
vcprank <- up_varvcp[order(up_varvcp$gene_symbol),]
petrank <- up_varpet[order(up_varpet$gene_symbol),]
ravrank <- up_varrav[order(up_varrav$gene_symbol),]

dev.new()
plot(up_varC9$variance....c9var, xlim = c(0,400), ylim = c(0,3000000), xlab = "Genes", ylab = "variance")
points(up_varftld$variance....ftldvar, col = "blue")
points(up_varvcp$variance....vcpvar, col = "green")
points(up_varsals$variance....salsvar, col = "orange")
points(up_varpet$variance....petvar, col = "red")
points(up_varrav$variance....ravvar, col = "magenta")



All_var <- data.frame(C9rank$gene_symbol,
                      C9orf72 = C9rank$rank,
                      sALS = salsrank$rank,
                      FTLD = ftldrank$rank,
                      VCP = vcprank$rank,
                      PET = petrank$rank,
                      RAV = ravrank$rank)

rownames(All_var) <- All_var$C9rank.gene_symbol
All_var[,1] <- NULL
All_var_rank <- All_var
All_var_rank$meanRank <- apply(All_var_sumrank, 1, mean)

RankVar <- apply(All_var, 1, var)
RankVar <- data.frame(row.names = rownames(All_var), 
                      variance <- RankVar)


setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
write.csv(RankVar, "upRankVar.csv")

# sumvar <- data.frame(row.names = rownames(All_var), 
#                      sumRank = All_var$sumRank)

nodetable <- read.csv("upGeneNode.csv")
nodetable <- as.data.frame(nodetable)
nodetable2 <- merge(nodetable, RankVar, by.x = "gene.name", by.y = 0)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
write.csv(nodetable2, "upGeneNode.csv")


############################################################################
###########DOWN VARIANCE#################################################
############################################################################


### PATIENT VARIANCE ####
options(scipen=999)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
rownames(C9) <- C9$Gene.Symbol
C9 <- C9[,12:19]

sals <- read.csv("sals_unique.csv")
rownames(sals) <- sals$Gene.Symbol
sals <- sals[,12:18]

ftld <- read.csv("ftld_unique.csv")
rownames(ftld) <- ftld$Gene.Symbol
ftld <- ftld[,17:32]

vcp <- read.csv("vcp_unique.csv")
rownames(vcp) <- vcp$Gene.Symbol
vcp <- vcp[,12:18]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
pet <- read.csv("PET_results_keepfiltering.csv")
rownames(pet) <- pet$hgnc_symbol
pet <- pet[,19:35]

rav <- read.csv("RAV_results_keepfiltering.csv")
rownames(rav) <- rav$hgnc_symbol
rav <- rav[,18:30]

## Calculate rowise variance
c9var <- apply(C9, 1, var)
salsvar <- apply(sals, 1, var)
ftldvar <- apply(ftld, 1, var)
vcpvar <- apply(vcp, 1, var)
petvar <- apply(pet, 1, var)
ravvar <- apply(rav, 1, var)

varC9 <- data.frame(row.names = rownames(C9), 
                    variance <- c9var)
varsals <- data.frame(row.names = rownames(sals), 
                      variance <- salsvar)
varftld <- data.frame(row.names = rownames(ftld), 
                      variance <- ftldvar)
varvcp <- data.frame(row.names = rownames(vcp), 
                     variance <- vcpvar)
varpet <- data.frame(row.names = rownames(pet), 
                     variance <- petvar)
varrav <- data.frame(row.names = rownames(rav), 
                     variance <- ravvar)

names <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/intersect_down_1.txt")
names <- names$V1 

up_varC9 <- subset(varC9, rownames(varC9) %in% names)
up_varsals<- subset(varsals, rownames(varsals) %in% names)
up_varftld <- subset(varftld, rownames(varftld) %in% names)
up_varvcp <- subset(varvcp, rownames(varvcp) %in% names)
up_varpet <- subset(varpet, rownames(varpet) %in% names)
up_varrav <- subset(varrav, rownames(varrav) %in% names)

up_varC9$gene_symbol <- rownames(up_varC9)
up_varC9 <- up_varC9[order(up_varC9$variance....c9var),]
up_varC9$rank <- 1:nrow(up_varC9)

up_varsals$gene_symbol <- rownames(up_varsals)
up_varsals <- up_varsals[order(up_varsals$variance....salsvar),]
up_varsals$rank <- 1:nrow(up_varsals)

up_varftld$gene_symbol <- rownames(up_varftld)
up_varftld <- up_varftld[order(up_varftld$variance....ftldvar),]
up_varftld$rank <- 1:nrow(up_varftld)

up_varvcp$gene_symbol <- rownames(up_varvcp)
up_varvcp <- up_varvcp[order(up_varvcp$variance....vcpvar),]
up_varvcp$rank <- 1:nrow(up_varvcp)

up_varpet$gene_symbol <- rownames(up_varpet)
up_varpet <- up_varpet[order(up_varpet$variance....petvar),]
up_varpet$rank <- 1:nrow(up_varpet)

up_varrav$gene_symbol <- rownames(up_varrav)
up_varrav <- up_varrav[order(up_varrav$variance....ravvar),]
up_varrav$rank <- 1:nrow(up_varrav)

C9rank <- up_varC9[order(up_varC9$gene_symbol),]
salsrank <- up_varsals[order(up_varsals$gene_symbol),]
ftldrank <- up_varftld[order(up_varftld$gene_symbol),]
vcprank <- up_varvcp[order(up_varvcp$gene_symbol),]
petrank <- up_varpet[order(up_varpet$gene_symbol),]
ravrank <- up_varrav[order(up_varrav$gene_symbol),]

dev.new()
plot(up_varC9$variance....c9var, xlim = c(0,400), ylim = c(0,3000000), xlab = "Genes", ylab = "variance")
points(up_varftld$variance....ftldvar, col = "blue")
points(up_varvcp$variance....vcpvar, col = "green")
points(up_varsals$variance....salsvar, col = "orange")
points(up_varpet$variance....petvar, col = "red")
points(up_varrav$variance....ravvar, col = "magenta")



All_var <- data.frame(C9rank$gene_symbol,
                      C9orf72 = C9rank$rank,
                      sALS = salsrank$rank,
                      FTLD = ftldrank$rank,
                      VCP = vcprank$rank,
                      PET = petrank$rank,
                      RAV = ravrank$rank)

rownames(All_var) <- All_var$C9rank.gene_symbol
All_var[,1] <- NULL
All_var_rank <- All_var
All_var_rank$meanRank <- apply(All_var_sumrank, 1, mean)

RankVar <- apply(All_var, 1, var)
RankVar <- data.frame(row.names = rownames(All_var), 
                    variance <- RankVar)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
write.csv(RankVar, "downRankVar.csv")

# sumvar <- data.frame(row.names = rownames(All_var), 
#                      sumRank = All_var$sumRank)

nodetable <- read.csv("downGeneNode.csv")
nodetable <- as.data.frame(nodetable)
nodetable2 <- merge(nodetable, RankVar, by.x = "gene.name", by.y = 0)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
write.csv(nodetable2, "downGeneNode.csv")
















# ########################
# ### VARIANCE ####
# options(scipen=999)
# 
# setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")
# 
# C9 <- read.csv("C9_unique.csv")
# rownames(C9) <- C9$Gene.Symbol
# C9 <- C9[,12:19]
# 
# sals <- read.csv("sals_unique.csv")
# rownames(sals) <- sals$Gene.Symbol
# sals <- sals[,12:18]
# 
# ftld <- read.csv("ftld_unique.csv")
# rownames(ftld) <- ftld$Gene.Symbol
# ftld <- ftld[,17:32]
# 
# vcp <- read.csv("vcp_unique.csv")
# rownames(vcp) <- vcp$Gene.Symbol
# vcp <- vcp[,12:18]
# 
# setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
# pet <- read.csv("PET_results_keepfiltering.csv")
# rownames(pet) <- pet$hgnc_symbol
# pet <- pet[,19:35]
# 
# rav <- read.csv("RAV_results_keepfiltering.csv")
# rownames(rav) <- rav$hgnc_symbol
# rav <- rav[,18:30]
# 
# ## Calculate rowise variance
# c9var <- apply(C9, 1, var)
# salsvar <- apply(sals, 1, var)
# ftldvar <- apply(ftld, 1, var)
# vcpvar <- apply(vcp, 1, var)
# petvar <- apply(pet, 1, var)
# ravvar <- apply(rav, 1, var)
# 
# varC9 <- data.frame(row.names = rownames(C9), 
#                     variance <- c9var)
# varsals <- data.frame(row.names = rownames(sals), 
#                       variance <- salsvar)
# varftld <- data.frame(row.names = rownames(ftld), 
#                       variance <- ftldvar)
# varvcp <- data.frame(row.names = rownames(vcp), 
#                      variance <- vcpvar)
# varpet <- data.frame(row.names = rownames(pet), 
#                      variance <- petvar)
# varrav <- data.frame(row.names = rownames(rav), 
#                      variance <- ravvar)
# 
# varC9$gene_symbol <- rownames(varC9)
# varC9 <- varC9[order(varC9$variance....c9var),]
# varC9$rank <- 1:nrow(varC9)
# 
# varsals$gene_symbol <- rownames(varsals)
# varsals <- varsals[order(varsals$variance....salsvar),]
# varsals$rank <- 1:nrow(varsals)
# 
# varftld$gene_symbol <- rownames(varftld)
# varftld <- varftld[order(varftld$variance....ftldvar),]
# varftld$rank <- 1:nrow(varftld)
# 
# varvcp$gene_symbol <- rownames(varvcp)
# varvcp <- varvcp[order(varvcp$variance....vcpvar),]
# varvcp$rank <- 1:nrow(varvcp)
# 
# varpet$gene_symbol <- rownames(varpet)
# varpet <- varpet[order(varpet$variance....petvar),]
# varpet$rank <- 1:nrow(varpet)
# 
# varrav$gene_symbol <- rownames(varrav)
# varrav <- varrav[order(varrav$variance....ravvar),]
# varrav$rank <- 1:nrow(varrav)
# 
# names <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/intersect_up_1.txt")
# names <- names$V1 
# up_varC9 <- subset(varC9, varC9$gene_symbol %in% names)
# up_varsals<- subset(varsals, varsals$gene_symbol %in% names)
# up_varftld <- subset(varftld, varftld$gene_symbol %in% names)
# up_varvcp <- subset(varvcp, varvcp$gene_symbol %in% names)
# up_varpet <- subset(varpet,varpet$gene_symbol %in% names)
# up_varrav <- subset(varrav, varrav$gene_symbol %in% names)
# 
# 
# C9rank <- up_varC9[order(up_varC9$gene_symbol),]
# salsrank <- up_varsals[order(up_varsals$gene_symbol),]
# ftldrank <- up_varftld[order(up_varftld$gene_symbol),]
# vcprank <- up_varvcp[order(up_varvcp$gene_symbol),]
# petrank <- up_varpet[order(up_varpet$gene_symbol),]
# ravrank <- up_varrav[order(up_varrav$gene_symbol),]
# 
# dev.new()
# plot(up_varC9$variance....c9var, xlim = c(0,400), ylim = c(0,3000000), xlab = "Genes", ylab = "variance")
# points(up_varftld$variance....ftldvar, col = "blue")
# points(up_varvcp$variance....vcpvar, col = "green")
# points(up_varsals$variance....salsvar, col = "orange")
# points(up_varpet$variance....petvar, col = "red")
# points(up_varrav$variance....ravvar, col = "magenta")
# 
# 
# 
# All_var <- data.frame(C9rank$gene_symbol,
#                       C9orf72 = C9rank$rank,
#                       sALS = salsrank$rank,
#                       FTLD = ftldrank$rank,
#                       VCP = vcprank$rank,
#                       PET = petrank$rank,
#                       RAV = ravrank$rank)
# 
# rownames(All_var) <- All_var$C9rank.gene_symbol
# All_var[,1] <- NULL
# 
# All_var$sumRank <- apply(All_var, 1, sum)
