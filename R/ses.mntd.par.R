#' Standardized effect size of MNTD
#'
#' Parallel calculation of standardized effect size of mean nearest taxon distances in communities. When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
#' 
#' Faster than ses.mntd from \code{picante} when there are many samples and taxa.
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
#' @keywords ses MNTD NTI
#' @return
#' A data frame of results for each community
#' \itemize{
#' \item{ntaxa}{Number of taxa in community}
#' \item{mntd.obs}{Observed mntd in community}
#' \item{mntd.rand.mean}{Mean mntd in null communities}
#' \item{mntd.rand.sd}{Standard deviation of mntd in null communities}
#' \item{mntd.obs.rank}{Rank of observed mntd vs. null communities}
#' \item{mntd.obs.z}{Standardized effect size of mntd vs. null communities (= (mntd.obs - mntd.rand.mean) / mntd.rand.sd, equivalent to -NRI)}
#' \item{mntd.obs.p}{P-value (quantile) of observed mntd vs. null communities (= mntd.obs.rank / runs + 1)}
#' \item{runs}{Number of randomizations}  
#' }
#' @import picante
#' @import doSNOW
#' @export

ses.mntd.par <- function(samp, dis, null.model = c("taxa.labels", "richness", 
                                                  "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
                                                  "trialswap"), abundance.weighted = FALSE, runs = 999, iterations = 1000, cores = 1){
    dis <- as.matrix(dis)
    mntd.obs <- mntd(samp, dis, abundance.weighted = abundance.weighted)
    null.model <- match.arg(null.model)
    N <- nrow(samp)
    
    # MNTD function on single sample
    mntd.single <- function(samp, dis, abundance.weighted, i){
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            diag(sample.dis) <- NA
            if (abundance.weighted) {
                mntds <- apply(sample.dis, 2, min, na.rm = TRUE)
                sample.weights <- samp[i, sppInSample]
                mntd <- weighted.mean(mntds, sample.weights)
            }
            else {
                mntd <- mean(apply(sample.dis, 2, min, na.rm = TRUE))
            }
        }
        else {
            mntd <- NA
        }
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
    mntd.rand <- foreach(i = seq_len(N), .options.snow = opts, .combine = cbind, .packages = "picante") %dopar% {
        
        rand.sub <- switch(null.model, 
                           taxa.labels = replicate(runs, mntd.single(samp, taxaShuffle(dis), abundance.weighted = abundance.weighted, i)), 
                           richness = replicate(runs, mntd.single(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted, i)), 
                           frequency = replicate(runs, mntd.single(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted, i)), 
                           sample.pool = replicate(runs, mntd.single(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted, i)), 
                           phylogeny.pool = replicate(runs, mntd.single(randomizeMatrix(samp, null.model = "richness"), taxaShuffle(dis), abundance.weighted, i)), 
                           independentswap = replicate(runs, mntd.single(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted, i)), 
                           trialswap = replicate(runs, mntd.single(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted, i)))
        
        return(rand.sub)
        
    }
    
    mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    mntd.obs.z <- (mntd.obs - mntd.rand.mean)/mntd.rand.sd
    mntd.obs.rank <- apply(X = rbind(mntd.obs, mntd.rand), MARGIN = 2,   FUN = rank)[1, ]
    mntd.obs.rank <- ifelse(is.na(mntd.rand.mean), NA, mntd.obs.rank)
    data.frame(ntaxa = specnumber(samp), 
               mntd.obs, 
               mntd.rand.mean, 
               mntd.rand.sd, 
               mntd.obs.rank, 
               mntd.obs.z, 
               mntd.obs.p = mntd.obs.rank/(runs + 1), runs = runs, row.names = row.names(samp))
}
