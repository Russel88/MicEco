#' Standardized effect size of MPD
#'
#' Parallel calculation of standardized effect size of mean pairwise distances in communities. When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Relative Index (NRI).
#' 
#' Faster than ses.mpd from \code{picante} when there are many samples and taxa.
#' @param samp Community data matrix with samples as rows
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param null.model Null model to use (see Details section for description)
#' @param abundance.weighted Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE)
#' @param runs Number of randomizations
#' @param iterations Number of iterations to use for each randomization (for independent swap and trial null models)
#' @param cores Number of cores to use for parallel computing
#' @details Currently implemented null models (arguments to null.model):
#' \itemize{
#'  \item taxa.labels - Shuffle distance matrix labels (across all taxa included in distance matrix) 
#'  \item richness - Randomize community data matrix abundances within samples (maintains sample species richness) 
#'  \item frequency - Randomize community data matrix abundances within species (maintains species occurence frequency) 
#'  \item sample.pool - Randomize community data matrix by drawing species from pool of species occurring in at least one community (sample pool) with equal probability   
#'  \item phylogeny.pool- Randomize community data matrix by drawing species from pool of species occurring in the distance matrix (phylogeny pool) with equal probability   
#'  \item independentswap - Randomize community data matrix with the independent swap algorithm (Gotelli 2000) maintaining species occurrence frequency and sample species richness  
#'  \item trialswap - Randomize community data matrix with the trial-swap algorithm (Miklos & Podani 2004) maintaining species occurrence frequency and sample species richness  
#' }
#' @keywords ses MPD NRI
#' @return
#' A data frame of results for each community
#' \itemize{
#' \item{ntaxa}{Number of taxa in community}
#' \item{mpd.obs}{Observed mpd in community}
#' \item{mpd.rand.mean}{Mean mpd in null communities}
#' \item{mpd.rand.sd}{Standard deviation of mpd in null communities}
#' \item{mpd.obs.rank}{Rank of observed mpd vs. null communities}
#' \item{mpd.obs.z}{Standardized effect size of mpd vs. null communities (= (mpd.obs - mpd.rand.mean) / mpd.rand.sd, equivalent to -NRI)}
#' \item{mpd.obs.p}{P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank / runs + 1)}
#' \item{runs}{Number of randomizations}  
#' }
#' @import picante
#' @import doSNOW
#' @export

ses.mpd.par <- function(samp, dis, null.model = c("taxa.labels", "richness", 
                                                  "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
                                                  "trialswap"), abundance.weighted = FALSE, runs = 999, iterations = 1000, cores = 1){
    dis <- as.matrix(dis)
    mpd.obs <- mpd(samp, dis, abundance.weighted = abundance.weighted)
    null.model <- match.arg(null.model)
    N <- nrow(samp)
    
    # MPD function on single sample
    mpd.single <- function(samp, dis, abundance.weighted, i){
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            if (abundance.weighted) {
                sample.weights <- t(as.matrix(samp[i, sppInSample, 
                                                   drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                                                                                      drop = FALSE])
                mpd <- weighted.mean(sample.dis, sample.weights)
            }
            else {
                mpd <- mean(sample.dis[lower.tri(sample.dis)])
            }
        }
        else {
            mpd <- NA
        }
        mpd
    }
    
    # Progress bar
    pb <- txtProgressBar(max = N, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Start parallel
    if(cores == 1) {
        registerDoSEQ() 
    } else {
        cl <- parallel::makeCluster(cores)
        registerDoSNOW(cl)
        on.exit(stopCluster(cl))
    }
    
    # Run parallel
    i <- NULL
    mpd.rand <- foreach(i = seq_len(N), .options.snow = opts, .combine = cbind, .packages = "picante") %dopar% {
        
        rand.sub <- switch(null.model, 
                           taxa.labels = replicate(runs, mpd.single(samp, taxaShuffle(dis), abundance.weighted = abundance.weighted, i)), 
                           richness = replicate(runs, mpd.single(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted, i)), 
                           frequency = replicate(runs, mpd.single(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted, i)), 
                           sample.pool = replicate(runs, mpd.single(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted, i)), 
                           phylogeny.pool = replicate(runs, mpd.single(randomizeMatrix(samp, null.model = "richness"), taxaShuffle(dis), abundance.weighted, i)), 
                           independentswap = replicate(runs, mpd.single(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted, i)), 
                           trialswap = replicate(runs, mpd.single(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted, i)))
        
        return(rand.sub)
        
    }
    
    mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
    mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2,   FUN = rank)[1, ]
    mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)
    data.frame(ntaxa = specnumber(samp), 
               mpd.obs, 
               mpd.rand.mean, 
               mpd.rand.sd, 
               mpd.obs.rank, 
               mpd.obs.z, 
               mpd.obs.p = mpd.obs.rank/(runs + 1), runs = runs, row.names = row.names(samp))
}
