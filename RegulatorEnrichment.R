library(biomaRt)
library (pathprint)
data(list = c("chipframe", "genesets","pathprint.Hs.gs","platform.thresholds", "pluripotents.frame"))
options(stringsAsFactors = FALSE)

pathways <- as.data.frame(pathprint.Hs.gs)

x <- grep("*(Netpath)", names(pathprint.Hs.gs))
netpath <- pathprint.Hs.gs[x]


nug <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NuggetGenes.txt")
GE <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/GeneXplain/NuggetUpstreamRegulators_genes.txt")
IPA <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/GeneXplain/IPAupstreamregulators.txt")

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
nug_entrez <- getBM(attributes =c("entrezgene", "hgnc_symbol"), filters="hgnc_symbol", values=nug,  mart=mart)
GE_entrez <- getBM(attributes =c("entrezgene", "hgnc_symbol"), filters="hgnc_symbol", values=GE,  mart=mart)
IPA_entrez <- getBM(attributes =c("entrezgene", "hgnc_symbol"), filters="hgnc_symbol", values=IPA,  mart=mart)

#Load file with all genes
library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

sym.genes <- t(sym.genes)

allgenes <- sym.genes[!duplicated(sym.genes),]

#####
hyperPathway <-
  function(genelist, geneset, Nchip) # produces a list of pathways enrichments calculated using the hypergeometric distribution
  {
    # Calculate p-values
    
    Nsig <- length(genelist)
    hyper<-as.data.frame(matrix(nrow = length(geneset), ncol = 1))
    rownames(hyper) <- names(geneset)
    colnames(hyper) <- c("p-value")	
    # determine p-value using the hypergeometric distribution, setting p = 1 if overlap = 0
    
    for (i in 1:length(geneset)){
      if (length(intersect(genelist, unlist(geneset[i]))) < 1){
        hyper[i,1]<-1}
      
      else if (length(intersect(genelist, unlist(geneset[i]))) > 0){
        hyper[i,1]<-phyper(length(intersect(genelist, unlist(geneset[i]))), Nsig, Nchip-Nsig, length(unlist(geneset[i])), lower.tail = FALSE)
      }
    }
    
    # adjust for multiple testing using Benjamini & Hochberg (fdr) correction
    hyper[,2]<-p.adjust(hyper[,1], method = "BH")		
    # Obtain list genes for each pathway
    
    overlap <- vector("list", 0)
    for (i in 1:length(geneset)){
      temp.overlap <- list(intersect(genelist, unlist(geneset[[i]])))
      overlap <- append(overlap, temp.overlap)			}
    
    names(overlap)<-rownames(hyper)
    
    # Count number of list genes and total genes in each pathway
    
    for (i in 1:length(geneset)){	
      hyper[i,3] <- length(overlap[[i]])
      hyper[i,4] <- length(geneset[[i]])
    }
    
    hyper[,5]<-rownames(hyper)
    hyper<-cbind((1:length(hyper[,1])), hyper)
    colnames(hyper)<-c("ID", "P-value", "BHadjP-value", "nGenes", "nPathway", "Name")						
    return(hyper)
  } # end of hyperPathway

# run script
pathwayEnrichment_IPA <- hyperPathway(
  genelist = IPA_entrez$entrezgene,
  geneset = netpath,
  Nchip = length(allgenes)
)

GEsig <- subset(pathwayEnrichment_GE, pathwayEnrichment_GE$`BHadjP-value` < 0.05)
IPAsig <- subset(pathwayEnrichment_IPA, pathwayEnrichment_IPA$`BHadjP-value` < 0.05)

sig <- intersect(rownames(GEsig), rownames(IPAsig))
