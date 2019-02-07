#' Inter-community mean nearest taxon distance
#'
#' Parallel calculation of MNTD (mean nearest taxon distance) separating taxa in two communities, a measure of phylogenetic beta diversity
#' @param comm Community data matrix with samples as rows
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param abundance.weighted Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE)
#' @param exclude.conspecifics Should conspecific taxa in different communities be exclude from MNTD calculations? (default = FALSE)
#' @param cores How many cores should be used for parallel computing'
#' @param progress Show a progress bar
#' @keywords MNTD bMNTD bNTI betaMNTD betaNTI
#' @return A distance matrix of MNTD values
#' @import picante
#' @import doSNOW
#' @export

comdistnt.par <- function (comm, dis, abundance.weighted = FALSE, exclude.conspecifics = FALSE, cores = 1, progress = TRUE) {
  
  dat <- match.comm.dist(comm, dis)
  comm <- dat$comm
  dis <- dat$dist
  N <- dim(comm)[1]
  comm <- decostand(comm, method = "total", MARGIN = 1)
  
  sppInSamples <- apply(comm,1,function(x) names(which(x > 0)))
  
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
  
  i <- NULL
  comdisnt <- foreach (i = 1:(N-1),.combine = rbind, .options.snow = opts) %dopar% {
    
    comdisnt.sub <- as.numeric(rep(NA,N))
    
    for (j in (i + 1):N) {
      
      sppInSample1 <- sppInSamples[[i]]
      sppInSample2 <- sppInSamples[[j]]
      
      if ((length(sppInSample1) >= 1) && (length(sppInSample2) >= 1)) {
        sample.dis <- dis[sppInSample1, sppInSample2, drop = FALSE]
        if (exclude.conspecifics) {
          sample.dis[sample.dis == 0] <- NA
        }
        sample1NT <- apply(sample.dis, 1, min, na.rm = TRUE)
        sample1NT[sample1NT == Inf] <- NA
        sample2NT <- apply(sample.dis, 2, min, na.rm = TRUE)
        sample2NT[sample2NT == Inf] <- NA
        
        if (abundance.weighted) {
          sample1.weights <- as.numeric(comm[i, sppInSample1])
          sample2.weights <- as.numeric(comm[j, sppInSample2])
          if (any(is.na(sample1NT))) {
            miss <- which(is.na(sample1NT))
            sample1NT <- sample1NT[-miss]
            sample1.weights <- sample1.weights[-miss]
            sample1.weights <- sample1.weights/sum(sample1.weights)
          }
          if (any(is.na(sample2NT))) {
            miss <- which(is.na(sample2NT))
            sample2NT <- sample2NT[-miss]
            sample2.weights <- sample2.weights[-miss]
            sample2.weights <- sample2.weights/sum(sample2.weights)
          }
          sampleNT <- c(sample1NT, sample2NT)
          sample.weights <- c(sample1.weights, sample2.weights)
          comdisnt.sub[j] <- weighted.mean(sampleNT, sample.weights,  na.rm = TRUE)
        }
        else {
          comdisnt.sub[j] <- mean(c(sample1NT, sample2NT), na.rm = TRUE)
        }
      }
      else {
        comdisnt.sub[j] <- NA
      }
    }
    return(comdisnt.sub)
  }
  if(cores != 1) stopCluster(cl)
  comdisnt <- rbind(comdisnt,rep(NA,N))
  
  rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
  return(as.dist(t(comdisnt)))
}
