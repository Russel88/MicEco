#' CLR transformation of community matrix, with multiplicative zero replacement
#'
#' This function performs a CLR transformation on the input.
#' Prior to transformation it does a multiplicative zero replacement, which adds a pseudocount, but corrects the non-zero abundances such that log-ratios between non-zero elements are unchanged after correction.
#' @param mat A community matrix with features as rows
#' @param delta The pseudocount to exchange zeroes with. Zero-correction is multiplicative such that the log-ratios between any entirely non-zero features will not be affected by the pseudocount. 
#' @keywords clr pseudocount
#' @return CLR transformed community matrix
#' @export

clr <- function(mat, delta=1){
  
  # Zero correction
  mat_zc <- apply(mat, 2, function(y) sapply(y,function(x) ifelse(x==0,delta,(1-(sum(y==0)*delta)/sum(y))*x)))
  
  if(any(mat_zc <= 0)) stop("Community matrix should only contain positive values")
  
  # CLR transformation
  gm_mean = function(x){
    if(any(x < 0, na.rm = TRUE)){
      stop("Negative values not allowed")
    }
    exp(mean(log(x)))
  }
  
  mat_gm <- apply(mat_zc,2,function(y) gm_mean(y))
  mat_log <- t(log(t(mat_zc)/mat_gm))

  rownames(mat_log) <- rownames(mat)
  colnames(mat_log) <- colnames(mat)

  return(mat_log)
}
