#' Run \code{UniFrac} multiple times in parallel and take the average
#'
#' With unrooted phylogenies \code{UniFrac} sets the root randomly on the tree. The position of the root affects the results.
#' The function runs UniFrac multiple times, with different roots, and takes the average to smooth potential bias.
#' @param physeq Phyloseq object. Required
#' @param R Number of times to repeat calculation. Default 100
#' @param seed Random seed for reproducibility. Default 42
#' @param cores Number of cores to use for parallel computing. Default 1 aka not parallel
#' @param ... Additional arguments passed to the \code{UniFrac} function
#' @return A distance object with the average \code{UniFrac} distances 
#' @import foreach parallel doSNOW phyloseq
#' @importFrom abind abind
#' @export
UniFrac.multi <- function(physeq, R = 100, seed = 42, cores = 1, ...){
  
  # Seeds for reproducibility
  seeds <- seed + 1:R
  
  # Parallel
  if(cores != 1){
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }
  
  # Progress bar
  pb <- txtProgressBar(max = R, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # UniFracs
  UFs <- foreach(i = 1:R, .packages = "phyloseq", .options.snow = opts) %dopar% {
    set.seed(seeds[i])
    UniFrac(physeq, ...)
  }

  # Average Unifracs
  UFs.mat <- lapply(UFs, as.matrix)  
  UF <- as.dist(apply(abind(UFs.mat, along = 3),c(1,2),mean))
  
}


