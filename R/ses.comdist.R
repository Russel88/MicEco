#' Standardized effect size of inter-community MPD (betaMPD, betaNRI)
#'
#' Standardized effect size of MPD (mean pairwise distance) separating taxa in two communities, a measure of phylogenetic beta diversity
#' @param samp Community data matrix
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param null.model Null model to use (see Details section for description)
#' @param abundance.weighted Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE)
#' @param runs Number of randomizations
#' @param iterations Number of iterations to use for each randomization (for independent swap and trial null models)
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
#' @keywords ses MPD bMPD bNRI betaMPD betaNRI
#' @return A list of results:
#' \itemize{
#'   \item ntaxa - Number of taxa in community
#'   \item comdist.obs - Observed mpd between community
#'   \item comdist.rand.mean - Mean mpd between null communities
#'   \item comdist.rand.sd - Standard deviation of mpd between null communities
#'   \item comdist.obs.rank - Rank of observed mpd vs. null mpd
#'   \item comdist.obs.z - Standardized effect size of mpd vs. null mpd (= (comdist.obs - comdist.rand.mean) / comdist.rand.sd, equivalent to -betaNRI)
#'   \item comdist.obs.p - P-value (quantile) of observed mpd vs. null communities (= comdist.obs.rank / runs + 1)
#'   \item runs - Number of randomizations
#' }
#' @import picante
#' @export

ses.comdist <- function (samp, dis, null.model = c("taxa.labels", "richness", 
                                                   "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
                                                   "trialswap"), abundance.weighted = FALSE, runs = 999, iterations = 1000){
  dis <- as.matrix(dis)
  comdist.obs <- as.matrix(comdist(samp, dis, abundance.weighted = abundance.weighted))
  null.model <- match.arg(null.model)
  comdist.rand <- switch(null.model, 
                         taxa.labels = replicate(runs, as.matrix(comdist(samp, taxaShuffle(dis), abundance.weighted)), simplify = FALSE), 
                         richness = replicate(runs, as.matrix(comdist(randomizeMatrix(samp, 
                                                                                      null.model = "richness"), dis, abundance.weighted)), simplify = FALSE), 
                         frequency = replicate(runs, as.matrix(comdist(randomizeMatrix(samp, 
                                                                             null.model = "frequency"), dis, abundance.weighted)), simplify = FALSE), 
                         sample.pool = replicate(runs, as.matrix(comdist(randomizeMatrix(samp, 
                                                                               null.model = "richness"), dis, abundance.weighted)), simplify = FALSE), 
                         phylogeny.pool = replicate(runs, as.matrix(comdist(randomizeMatrix(samp, 
                                                                                  null.model = "richness"), taxaShuffle(dis), abundance.weighted)), simplify = FALSE), 
                         independentswap = replicate(runs, as.matrix(comdist(randomizeMatrix(samp, 
                                                                                   null.model = "independentswap", iterations), dis, abundance.weighted)), simplify = FALSE), 
                         trialswap = replicate(runs,as.matrix(comdist(randomizeMatrix(samp, 
                                                                            null.model = "trialswap", iterations), dis, abundance.weighted)), simplify = FALSE))
  
  comdist.rand.mean <- apply(X = simplify2array(comdist.rand), MARGIN = 1:2, FUN = mean, na.rm = TRUE)
  
  comdist.rand.sd <- apply(X = simplify2array(comdist.rand), MARGIN = 1:2, FUN = sd, na.rm = TRUE)
  
  comdist.obs.z <- (comdist.obs - comdist.rand.mean)/comdist.rand.sd
  
  comdist.obs.rank <- apply(X = simplify2array(c(list(comdist.obs),comdist.rand)), MARGIN = 1:2, FUN = rank)[1,,]
  comdist.obs.rank <- ifelse(is.na(comdist.rand.mean), NA, comdist.obs.rank)
  diag(comdist.obs.rank) <- NA
  
  comdist.obs.p = comdist.obs.rank/(runs + 1)
  
  list(ntaxa = specnumber(samp), comdist.obs = comdist.obs, comdist.rand.mean = comdist.rand.mean, 
       comdist.rand.sd = comdist.rand.sd, comdist.obs.rank = comdist.obs.rank, comdist.obs.z = comdist.obs.z, comdist.obs.p = comdist.obs.p, runs = runs)
}

