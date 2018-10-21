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
  
  # Zero correction
  otu_zc <- apply(otu, 2, function(y) sapply(y,function(x) ifelse(x==0,delta,(1-(sum(y==0)*delta)/sum(y))*x)))
  
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
  mat <- matrix(0,nrow=nrow(otu),ncol=nrow(otu))
  for(i in 1:(nrow(otu)-1)){
    for(j in (i+1):nrow(otu)){
      mat[i,j] <- (2*cov(otu_log[i,],otu_log[j,]))/(var(otu_log[i,])+var(otu_log[j,]))
    }
  }
  rownames(mat) <- rownames(otu)
  colnames(mat) <- rownames(otu)
  mat <- t(mat)+mat
  diag(mat) <- 1
  return(mat)
}
