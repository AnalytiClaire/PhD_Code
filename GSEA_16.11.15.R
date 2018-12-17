
###START HERE #####

setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/C9orf72_LCM/") #set working directory
exp <- read.csv("eset_NineP_150612_exprs.csv") #import expression file
row.names(exp) <- exp[,1] #make Probe IDs row names if not already done
exp[,1] <- NULL #remove column containing probe IDs

library (WGCNA)
options(stringsAsFactors = FALSE)

DATA <- exp 
DATA <- t(DATA) #transpose data (this is required for Biomart)


####GSEA

### connect to the biomart database that is correct for your data ###
library (biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org") 
x <- colnames(DATA) #create vector containing probe IDs


### Retrive the specified attributes from the above database ###

mart_attribute <- listAttributes(mart)
# mart_filter <- listFilters(mart)



### Create array of attributes required for your data set ### MAKE SURE YOU GET THE INFO OF YOUR PROBE IDs
mart_back <- getBM(attributes=c("hgnc_symbol","affy_hg_u133_plus_2"), 
                   filters = "affy_hg_u133_plus_2", values = x, mart = mart)

e2 <- mart_back$hgnc_symbol #take gene label column
e3 <- e2[!is.na(e2)]
e3 <- which(!(e2 == "")) #remove any rows with blank cells
mart_back1 <- mart_back[c(e3),] #apply to data frame





#e4 <- mart_back1[,1] #not sure what this is for

# #e2 <- mart_back[,2]
# e2 <- mart_back[,1]
# e3 <- which (!is.na(e2))
# mart_back1 <- mart_back[c(e3),]
# #e4 <- mart_back1[,2]
# e4 <- mart_back1[,1]
# 
# c <- vector(length=length (mart_back1[,3]))
# c <- vector(length=length (mart_back1[,3]))





### Assign hgnc symbol as column names of expression data ###
#This looks at the ID probes from the expression file and links them to the corresponding 
#hgnc symbol

#This will take a decent amount of time, so don't worry if it runs for 30 minutes+

DATA1 <- t(DATA)

for (i in 1:length(mart_back1$hgnc_symbol)) 
{
  c1 <- which (x %in% mart_back1$hgnc_symbol)
  c1 <- c1[1]
  #colnames (DATA1)[c(c1)] <- mart_back1[i,2]
  colnames (DATA1)[c(c1)] <- mart_back1[i,1]
}

# DATA1 <- cbind(Row.Names = rownames(DATA1), DATA1)
# 
# DATA1 <- merge(mart_back1, DATA1, by.x = "affy_hg_u133_plus_2", by.y = "Row.Names" )
# DATA1[,1] <- NULL

#agg.data <- aggregate(DATA1,by=list(DATA1$hgnc_symbol),mean)

#rownames(DATA1) <- DATA1[,1]

t <- array (dim =c(length(DATA[,1]), length (DATA[1,]), d))


####I highly recommend saving the environment as you don't want to have to run this twice ####





# #### edit background for GSEA
# b <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA-P-R/GeneSetDatabases/c5.bp.v5.0.symbols.gmt", 
#               header=FALSE, fill=TRUE)
# ba <- matrix (nrow=nrow(b), ncol=2)
# bb <- matrix (nrow=nrow(b), ncol=2000)
# 
# 
# dx <- colnames(DATA1)
# dx1 <- grep("ENST", dx)
# dx2 <- dx [-c(dx1)]
# 
# 
# for (i in 1:nrow(b))
# {
# b1 <- b[i,1]
# b1 <- strsplit (b1, split="\t")
# b1 <- b1[[1]]
# ba[i,] <- b1[1:2]
# 
# b2 <- b1[3:length(b1)]                 
# b3 <- which (b2 %in% dx2)
# if (length(b3)>0) {b2 <- b2[c(b3)]}
# if (length(b3)==0) {b2 <- NA}
# bb[i,1:length(b2)] <- b2
# bb[i,(length(b2)+1):ncol(bb)] <- NA
# 
# }
# 
# gmt <- cbind(ba,bb)
# write.table (gmt, file= "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA-P-R/GeneSetDatabases/c5.bp.v5.0.symbols.calib.gmt", 
#              row.names=FALSE, sep = "\t", col.names=FALSE, append = F, na="", quote=FALSE)

#############################currently need to remove NAs in excel #############

#load  ("C:/Users/sann1458/Dropbox/Transcriptome studies/Win260914/290115.RData")





###create gct file###

t2 <- t(DATA1) #t2 needs samples as row names and gene IDs as column names

# t1 <- grep("ENST", colnames(t)) #This is for when you are using data with ensembl transcript IDs
# t2 <- t [,-c(t1)]

t1a <- which (colnames(t2) %in% "") #take column names in t2 that contain blanks
if (length (t1a)>0) {t2 <- t2 [,-c(t1a)]} #if there are any blanks, remove them (this will be 0 if you do not)

t2.colnames <- toupper(colnames(t2)) #convert all column names of t2 to uppercase

#create an empty matrix with same rows as t2, but with only the same number of columns as unique names
t2.agg <- matrix (nrow=length(t2[,1]), ncol=length(unique(t2.colnames))) 
row.names (t2.agg) <- row.names (t2) #make the row names the same as t2
colnames (t2.agg) <- sort(unique(t2.colnames)) #make the colum names = to unique t2 column names and sort

#This aggregates the t2 columns by name, takes the mean, and then applies it back to the data frame 
for (i in 1:length (t2[,1]))
{
agg <- aggregate(t2[i,],by=list(t2.colnames),mean)
t2.agg[i,] <- agg[,2]
}



#j="All."
t.agg1 <- t(t2.agg)
t.agg2 <- cbind (row.names(t.agg1), rep(NA,length(t.agg1[,1])), t.agg1)
colnames(t.agg2)[1:2] <- c("NAME", "DESCRIPTION")
#t.agg2 <- data.matrix(t.agg2)
#mode(t.agg2) <- "numeric"
write("#1.2", file = paste("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/FTLD/", "FTLD.gct", sep=""), append = F)
write(paste(dim(t.agg2)[1], dim(t.agg1)[2], sep = "\t"),file = paste("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/FTLD/", "FTLD.gct", sep=""), append = T)
write.table (t.agg2, file= paste("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/FTLD/", "FTLD.gct", sep=""), row.names=FALSE, append = T, sep = "\t", col.names=TRUE, quote=FALSE)


#Create cls file
write(paste("24", "2", "1", sep = " "), file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/FTLD/FTLD_pheno.cls", append = F)
write(paste("#", "PAT", "CON", sep = " "),file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/FTLD/FTLD_pheno.cls", append = T)
write(paste("0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0",
            "1","1","1","1","1","1","1","1", sep = " "),file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA/FTLD/FTLD_pheno.cls", append = T)




###AT THIS POINT YOU CAN USE THE GSEA JAVA PLATFORM AS THE R SCRIPT IS NO LONGER SUPPORTED
                 
GSEA.program.location <- ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA-P-R/GSEA.1.0.R")   #  R source program (change pathname to the right location in local machine)
source(GSEA.program.location, verbose=T, max.deparse.length=9999)


#dir.create(paste("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/CHMP2B/GSEA_NO"))

GSEA(                                                                    # Input/Output Files :-------------------------------------------
 input.ds =  paste(""), # Input gene expression Affy dataset file in RES or GCT format
 input.cls = "", # Input class vector (phenotype) file in CLS format
 gs.db =     "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/GSEA-P-R/GeneSetDatabases/c5.bp.v5.0.symbols.gmt", # Gene set database in GMT format
 output.directory      = paste(",""/", sep=""),        # Directory where to store output and results (default: "")
#  Program parameters :-------------------------------------------------------------------------------------------------------------------------
 doc.string            = "",   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
 non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
 reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
 nperm                 = 1000,           # Number of random permutations (default: 1000)
 weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
 nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
 fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
 fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
 topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
 adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
 gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
 gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
 reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
 preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
 random.seed           = 3338,            # Random number generator seed. (default: 123456)
 perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
 fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
 replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
 save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
 OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
 use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
)
#-----------------------------------------------------------------------------------------------------------------------------------------------

results0 <- matrix(unlist(read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/FUS_SALS_LCM_CELfiles/sALS_GSEA_Output/sALS_GSEA.SUMMARY.RESULTS.REPORT.0.txt", sep="",fill=T)),
                   ncol=max(count.fields("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/FUS_SALS_LCM_CELfiles/sALS_GSEA_Output/sALS_GSEA.SUMMARY.RESULTS.REPORT.0.txt", sep="")),byrow=F)
results1 <- matrix(unlist(read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/FUS_SALS_LCM_CELfiles/sALS_GSEA_Output/sALS_GSEA.SUMMARY.RESULTS.REPORT.1.txt", sep="",fill=T)),
                   ncol=max(count.fields("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/FUS_SALS_LCM_CELfiles/sALS_GSEA_Output/sALS_GSEA.SUMMARY.RESULTS.REPORT.1.txt", sep="")),byrow=F)

print (results0[1:10,])
print (results1[1:10,])











