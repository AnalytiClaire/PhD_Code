function (genelist, geneset, Nchip, ogGeneset) 
{
  Nsig <- length(genelist)
  hyper <- as.data.frame(matrix(nrow = length(geneset), ncol = 1))
  rownames(hyper) <- names(geneset)
  colnames(hyper) <- c("p-value")
  for (i in 1:length(geneset)) {
    if (length(intersect(genelist, unlist(geneset[i]))) < 
        1) {
      hyper[i, 1] <- 1
    }
    else if (length(intersect(genelist, unlist(geneset[i]))) > 
             0) {
      hyper[i, 1] <- phyper(length(intersect(genelist, 
                                             unlist(geneset[i]))), Nsig, Nchip - Nsig, length(unlist(ogGeneset[i])), 
                            lower.tail = FALSE)
    }
  }
  hyper[, 2] <- p.adjust(hyper[, 1], method = "BH")
  overlap <- vector("list", 0)
  for (i in 1:length(geneset)) {
    temp.overlap <- list(intersect(genelist, unlist(geneset[[i]])))
    overlap <- append(overlap, temp.overlap)
  }
  names(overlap) <- rownames(hyper)
  for (i in 1:length(geneset)) {
    hyper[i, 3] <- length(overlap[[i]])
    hyper[i, 4] <- length(geneset[[i]])
  }
  hyper[, 5] <- rownames(hyper)
  hyper <- cbind((1:length(hyper[, 1])), hyper)
  colnames(hyper) <- c("ID", "P-value", "BHadjP-value", "nGenes", 
                       "nPathway", "Name")
  return(hyper)
}