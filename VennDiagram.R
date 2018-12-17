setwd(dir = "/Users/clairegreen/Documents/PhD/LRRK2_PARKIN_Project/AltAnalyzeDEGs")

library(VennDiagram)

a <- read.table(file = "man_VS_control (adjp_0.05).txt", header = FALSE)
a <- a$V1
b <- read.table(file = "man_vs_nonman (rawp_0.001).txt", header = FALSE)
b <- b$V1
c <- read.table(file = "Non-man_VS_control (adjp_0.05).txt", header = FALSE)
c <- c$V1

a <- unique(a)
b <- unique(b)
c <- unique(c)





results <- read.csv("LRRK2_manvsnonman.csv", na.strings = c("", "NA)"))
results <- as.list(results)
# results <- lapply(results, function(x) x[!duplicated(x)])
results <- lapply(results, function(x) x[!is.na(x)])


overlap <- calculate.overlap(results)

grid.newpage()
draw.triple.venn(area1 = 510, 
                 area2 = 1308, 
                 area3 = 3400, 
                 n12 = 101, 
                 n23 = 293, 
                 n13 = 294, 
                 n123 = 61, 
                 euler.d =TRUE, scaled = TRUE,
                 category = c("CG (DESeq2)", "SA (DESeq2)", "LC (BitSeq"), lty = "blank",
                 fill = c("magenta3", "dodgerblue3", "yellow3"), cat.dist = c(-0.06, -0.055, -0.04), cex = 1.5, cat.cex = c(2,2,2))


MVC_NMVC <- intersect(results$man.vs.con,results$nonman.vs.con)
MVC_MVNM <- intersect(results$man.vs.con,results$man.vs.nonman)
NMVC_MVNM <- intersect(results$nonman.vs.con,results$man.vs.nonman)

MVC_unique <- setdiff(MVC_NMVC, Mmc)
MVNM_unique <- setdiff(results$man.vs.con,results$man.vs.nonman)
NMVC_unique <- setdiff(results$nonman.vs.con,results$man.vs.nonman)

write.table(MVC_MVNM, file = "overlap_MVC_NVNM.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(MVC_NMVC, file = "overlap_MVC_NMVC.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(NMVC_MVNM, file = "overlap_NMVC_MVNM.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)


venn.plot <- draw.triple.venn(20, 40, 60, 0, 0, 0, 0,
                              c("First", "Second", "Third"), sep.dist = 0.1, rotation.degree = 30)

