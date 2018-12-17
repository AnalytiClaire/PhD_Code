# library(devtools)
# install_github("emreg00/pepper")

# Data set specific parameters
geo.id <- "GSE68605" # GEO id of the data set
probe.conversion <- "ENTREZ_GENE_ID" # column name for gene id mapping
conversion.map <- NULL # probe to gene mapping, if NULL uses the mapping in the data set
conversion.mapping.function <- NULL # modify probe names using this function 
sample.mapping.column <- "characteristics_ch1" # column to use for sample mapping
geo.id.sub = NULL # the platform to use if there are multiple platform annotations
reprocess <- "affy" # reprocessing type for raw data
output.dir <- "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/PePPer"

# Get the expression and sample mapping info from the reprocessed data set
d <- fetch.expression.data(geo.id, 
                           sample.mapping.column = sample.mapping.column, 
                           do.log2 = F, 
                           probe.conversion = probe.conversion)


expr <- d$expr
sample.mapping <- d$sample.mapping

#GET GROUPWISE DE
# states.case <- c("Patient")
# states.control <- c("Control") 
# sample.mapping <- convert.sample.mapping.to.case.control(sample.mapping, states.case, states.control)
sample.mapping$type = c(c(rep("case", 8)), c(rep("control",3)))

adjust.method <- 'BH'
fdr.cutoff <- 0.05
out.file <- "de.dat"
de <- find.de.genes(expr, sample.mapping, c("case", "control"), method="limma", adjust.method=adjust.method, cutoff=fdr.cutoff, functional.enrichment="kegg") 
de <- de[abs(de$logFC)>=1,]

#GET INDIVIDUAL DE
# Get z scores
out.file <- "z.dat"
cutoff <- 2.5
z = get.z.matrix(expr, sample.mapping, method="mean", out.file=out.file)
indices <- apply(abs(z), 2, function(x) { which(x >= cutoff)})
geneids <- lapply(indices, function(x) { rownames(z)[x] })
