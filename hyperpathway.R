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