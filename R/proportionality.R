#' Calculate proportionality between all pairs of OTUs, for network analysis
#'
#' Calculate proportionality as proposed by Lovell et al. 2016 Proportionality: a valid alternative to correlation for relative data
#' @param x An otu-table with OTUs as rows OR a phyloseq object
#' @param delta The pseudocount to exchange zeroes with. Zero-correction is multiplicative such that the proportionality between any entirely non-zero OTUs will not be affected by the pseudocount. 
#' @keywords network association
#' @return The network association matrix with proportionality values
#' @import phyloseq
#' @export

proportionality <- function(x,delta=1){
  
  if(class(x)!="phyloseq") otu <- x else {
    if(taxa_are_rows(x)) otu <- otu_table(x)
    if(!taxa_are_rows(x)) otu <- t(otu_table(x)) 
  }
  
  mat <- matrix(nrow=nrow(otu),ncol=nrow(otu)) 
  
  # Zero correction
  otu_zc <- matrix(nrow=nrow(otu),ncol=ncol(otu))  
  for(k in 1:ncol(otu)){
    samp <- as.data.frame(otu[,k])
    z <- sum(samp==0)
    col_new <- apply(samp,1,function(x) ifelse(x==0,delta,(1-(z*delta)/sum(samp))*x))
    otu_zc[,k] <- col_new
  }
  
  if(any(otu_zc <= 0)) stop("OTU table should only contain positive values")
  
  # CLR transformation
  gm_mean = function(x){
    if(any(x < 0, na.rm = TRUE)){
      stop("Negative values not allowed")
    }
    exp(mean(log(x)))
  }
  
  otu_gm <- apply(otu_zc,2,function(y) gm_mean(y))
  otu_log <- t(log(t(otu_zc)/otu_gm))
  
  # Proportionality
  for(i in 1:nrow(otu)){
    for(j in 1:nrow(otu)){
      if(i==j) prop <- 1
      if(i != j) prop <- (2*cov(otu_log[i,],otu_log[j,]))/(var(otu_log[i,])+var(otu_log[j,]))
      mat[i,j] <- prop
    }
  }
  rownames(mat) <- rownames(otu)
  colnames(mat) <- rownames(otu)
  return(mat)
}
