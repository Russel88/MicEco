#' Bootstrap of z-matrix 
#' 
#' Bootstrap of z-matrix from ses.comdist(2), ses.comdistnt and ses.UniFrac.
#' Estimate distribution for each group and combination of groups.
#' Also tests if distribution is different from a normal distribution with mean zero and sd equal to bootstrapped distribution.
#' @param zmat Symmetric matrix with z-values
#' @param sampleGroups Vector with the grouping of the samples in the same order as in zmat
#' @param R Number of bootstraps
#' @param probs The quantiles to return. Default is 95% quantiles and the median.
#' @return Quantiles given by probs, and adjusted p-values testing whether estimate is above or below zero
#' @export

ses.bootstrap <- function(zmat, sampleGroups, R = 10000, probs = c(0.025, 0.5, 0.975)){
  
  ### Within
  # Subset matrices
  mats.w <- lapply(unique(sampleGroups), function(x) zmat[sampleGroups == x,sampleGroups == x])

  # NA upper tris
  for(i in 1:length(mats.w)){
    mats.w[[i]][upper.tri(mats.w[[i]], diag = TRUE)] <- NA
  }
  rm(i)
  
  # Sample rows
  samp_row <- lapply(1:length(mats.w), function(y) lapply(1:R,function(x) sample(1:nrow(mats.w[[y]]), replace = TRUE)))
  samp_col <- lapply(1:length(mats.w), function(y) lapply(1:R,function(x) sample(1:ncol(mats.w[[y]]), replace = TRUE)))
  
  # Bootstrap
  boot_within <- lapply(1:length(mats.w), function(x) sapply(1:R, function(y) mean(mats.w[[x]][samp_row[[x]][[y]],samp_col[[x]][[y]]],na.rm=TRUE)))
  
  
  ### Between
  combis <- combn(unique(sampleGroups),2)
  
  # Subset matrices
  mats.b <- lapply(1:ncol(combis), function(x) zmat[sampleGroups == combis[1,x],sampleGroups == combis[2,x]])
  
  # Sample rows
  samp_row.b <- lapply(1:length(mats.b), function(y) lapply(1:R,function(x) sample(1:nrow(mats.b[[y]]), replace = TRUE)))
  samp_col.b <- lapply(1:length(mats.b), function(y) lapply(1:R,function(x) sample(1:ncol(mats.b[[y]]), replace = TRUE)))
  
  # Bootstrap
  boot_between <- lapply(1:length(mats.b), function(x) sapply(1:R, function(y) mean(mats.b[[x]][samp_row.b[[x]][[y]],samp_col.b[[x]][[y]]])))
  
  
  # Combine
  boots <- c(boot_within,boot_between)
  
  # Output
  CIs <- t(sapply(boots, function(x) quantile(x, probs = probs)))
  rownames(CIs) <- c(unique(sampleGroups),apply(combis,2,function(x) paste(x, collapse = "_")))
  
  pvals.below <- sapply(boots, function(x) (sum(x >= rnorm(length(x),mean=0,sd=sd(x)))+1)/(R+1))
  pvals.above <- sapply(boots, function(x) (sum(x <= rnorm(length(x),mean=0,sd=sd(x)))+1)/(R+1))
  
  CIs <- as.data.frame(CIs)
  
  CIs$pval.adj.below <- p.adjust(pvals.below, method = "fdr")
  CIs$pval.adj.above <- p.adjust(pvals.above, method = "fdr")
  
  return(CIs)
}
