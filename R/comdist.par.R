#' Inter-community mean pairwise distance
#'
#' Parallel calculation of MPD (mean pairwise distance) separating taxa in two communities, a measure of phylogenetic beta diversity
#' @param comm Community data matrix with samples as rows
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param abundance.weighted Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE)
#' @param cores How many cores should be used for parallel computing
#' @param progress Show a progress bar
#' @keywords MPD bMPD bNRI betaMPD betaNRI
#' @return A distance matrix of MPD values
#' @import picante
#' @import doSNOW
#' @export

comdist.par <- function (comm, dis, abundance.weighted = FALSE, cores = 1, progress = TRUE) {
  
  x <- as.matrix(comm)
  dat <- match.comm.dist(comm, dis)
  x <- dat$comm
  dis <- as.matrix(dat$dist)
  if (!abundance.weighted) {
    x <- decostand(x, method = "pa")
  }
  N <- dim(x)[1]
  S <- dim(x)[2]
  x <- decostand(x, method = "total", MARGIN = 1)
  
  if(progress){
    pb <- txtProgressBar(max = (N-1), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }
 
  if(cores == 1) {
    registerDoSEQ() } 
  else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  }
  
  l <- NULL 
  comdist <- foreach(l = 1:(N - 1),.combine = cbind, .options.snow = opts) %dopar% {
    comdist.sub <- as.numeric(rep(NA,N))
    for (k in 2:N) {
      comdist.sub[k] <- sum(dis * outer(as.vector(t(x[k,])), as.vector(t(x[l, ]))))
    }
    return(comdist.sub)
  }
  if(cores != 1) stopCluster(cl)
  comdist <- cbind(comdist,rep(NA,N))
  row.names(comdist) <- row.names(x)
  colnames(comdist) <- row.names(x)
  return(as.dist(comdist))
}
