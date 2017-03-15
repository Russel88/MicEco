#' Standardized effect size of inter-community MNTD (betaMNTD, betaNTI)
#'
#' Standardized effect size of MNTD (mean nearest taxon distance) separating taxa in two communities, a measure of phylogenetic beta diversity
#' @param samp Community data matrix with samples as rows
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param null.model Null model to use (see Details section for description)
#' @param abundance.weighted Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE)
#' @param exclude.conspecifics Should conspecific taxa in different communities be exclude from MNTD calculations? (default = FALSE)
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
#' @keywords ses MNTD bMNTD bNTI betaMNTD betaNTI
#' @return A list of results:
#' \itemize{
#'   \item ntaxa - Number of taxa in community
#'   \item comdistnt.obs - Observed mntd between communities
#'   \item comdistnt.rand.mean - Mean mntd between null communities
#'   \item comdistnt.rand.sd - Standard deviation of mntd between null communities
#'   \item comdistnt.obs.rank - Rank of observed mntd vs. null mntd
#'   \item comdistnt.obs.z - Standardized effect size of mntd vs. null mntd (= (comdistnt.obs - comdistnt.rand.mean) / comdistnt.rand.sd, equivalent to -betaNTI)
#'   \item comdistnt.obs.p - P-value (quantile) of observed mntd vs. null communities (= comdistnt.obs.rank / runs + 1)
#'   \item runs - Number of randomizations
#' }
#' @import picante
#' @export

ses.comdistnt <- function (samp, dis, null.model = c("taxa.labels", "richness", 
                                                     "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
                                                     "trialswap"), abundance.weighted = FALSE, exclude.conspecifics = FALSE, runs = 999, iterations = 1000){
  dis <- as.matrix(dis)
  comdistnt.obs <- as.matrix(comdistnt(samp, dis, abundance.weighted = abundance.weighted, exclude.conspecifics = exclude.conspecifics))
  null.model <- match.arg(null.model)
  comdistnt.rand <- switch(null.model, taxa.labels = replicate(runs, 
                                                                 as.matrix(comdistnt(samp, taxaShuffle(dis), abundance.weighted = abundance.weighted, exclude.conspecifics = exclude.conspecifics)), simplify = FALSE), 
                           richness = replicate(runs, as.matrix(comdistnt(randomizeMatrix(samp, 
                                                                                          null.model = "richness"), dis, abundance.weighted, exclude.conspecifics = exclude.conspecifics)), simplify = FALSE), 
                           frequency = replicate(runs, as.matrix(comdistnt(randomizeMatrix(samp, 
                                                                                 null.model = "frequency"), dis, abundance.weighted, exclude.conspecifics = exclude.conspecifics)), simplify = FALSE), 
                           sample.pool = replicate(runs, as.matrix(comdistnt(randomizeMatrix(samp, 
                                                                                   null.model = "richness"), dis, abundance.weighted, exclude.conspecifics = exclude.conspecifics)), simplify = FALSE), 
                           phylogeny.pool = replicate(runs, as.matrix(comdistnt(randomizeMatrix(samp, 
                                                                                      null.model = "richness"), taxaShuffle(dis), abundance.weighted, exclude.conspecifics = exclude.conspecifics)), simplify = FALSE), 
                           independentswap = replicate(runs, as.matrix(comdistnt(randomizeMatrix(samp, 
                                                                                       null.model = "independentswap", iterations), dis, abundance.weighted, exclude.conspecifics = exclude.conspecifics)), simplify = FALSE), 
                           trialswap = replicate(runs,as.matrix(comdistnt(randomizeMatrix(samp, 
                                                                                null.model = "trialswap", iterations), dis, abundance.weighted, exclude.conspecifics = exclude.conspecifics)), simplify = FALSE))
  
  comdistnt.rand.mean <- apply(X = simplify2array(comdistnt.rand), MARGIN = 1:2, FUN = mean, na.rm = TRUE)
  
  comdistnt.rand.sd <- apply(X = simplify2array(comdistnt.rand), MARGIN = 1:2, FUN = sd, na.rm = TRUE)
  
  comdistnt.obs.z <- (comdistnt.obs - comdistnt.rand.mean)/comdistnt.rand.sd
  
  comdistnt.obs.rank <- apply(X = simplify2array(c(list(comdistnt.obs),comdistnt.rand)), MARGIN = 1:2, FUN = rank)[1,,]
  comdistnt.obs.rank <- ifelse(is.na(comdistnt.rand.mean), NA, comdistnt.obs.rank)
  diag(comdistnt.obs.rank) <- NA
  
  comdistnt.obs.p <- comdistnt.obs.rank/(runs + 1)
  
  list(ntaxa = specnumber(samp), comdistnt.obs = comdistnt.obs, comdistnt.rand.mean = comdistnt.rand.mean, 
       comdistnt.rand.sd = comdistnt.rand.sd, comdistnt.obs.rank = comdistnt.obs.rank, comdistnt.obs.z = comdistnt.obs.z, comdistnt.obs.p = comdistnt.obs.p, runs = runs)
}
