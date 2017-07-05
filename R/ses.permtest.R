#' Permutation test of z-matrix 
#' 
#' Permutation test of z-matrix from ses.comdist(2), ses.comdistnt and ses.UniFrac.
#' Test if means within or between groups is higher or lower than the overall mean of the matrix
#' @param zmat Symmetric matrix with z-values
#' @param sampleGroups Vector with the grouping of the samples in the same order as in zmat
#' @param R Number of permutations
#' @return A dataframe with p-values (pval), fdr corrected p-values (pval.adj), and average z-values (avg) for each group and combination of groups.
#' @export

ses.permtest <- function(zmat, sampleGroups, R = 10000){
  
  # Permute labels
  perms <- lapply(1:R,function(x) sample(sampleGroups))
  
  ### Within groups
  # Real means
  real.mean <- lapply(unique(sampleGroups), function(x) mean(zmat[sampleGroups == x,sampleGroups == x][lower.tri(zmat[sampleGroups == x,sampleGroups == x])]))
  
  # Means for permuted groups
  perm.mean <- lapply(unique(sampleGroups), function(x) sapply(perms, function(y) mean(zmat[y == x,y == x][lower.tri(zmat[y == x,y == x])])))
  
  ### Between groups
  combis <- combn(unique(sampleGroups),2)
  
  # Real means
  real.mean.c <- lapply(1:ncol(combis), function(x) mean(zmat[sampleGroups == combis[1,x],sampleGroups == combis[2,x]]))
  
  # Means for permuted groups
  perm.mean.c <- lapply(1:ncol(combis), function(x) sapply(perms, function(y) mean(zmat[y == combis[1,x],y == combis[2,x]])))
  
  # Put together
  real.means <- c(real.mean,real.mean.c)
  perm.means <- c(perm.mean,perm.mean.c)
  
  # P-values
  pvals <- sapply(1:length(real.means), function(x) (sum(abs(perm.means[[x]]) >= abs(mean(real.means[[x]])))+1)/(R+1))
  
  # Output
  names(pvals) <- c(unique(sampleGroups),apply(combis,2,function(x) paste(x, collapse = "_")))
  df <- data.frame(pval = pvals,
                   pval.adj = p.adjust(pvals, method = "fdr"),
                   avg = as.numeric(real.means))
  
  return(df)
}