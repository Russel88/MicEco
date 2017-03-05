#' Standardized effect size of inter-community MPD (betaMPD, betaNRI)
#'
#' Standardized effect size of MPD (mean pairwise distance) separating taxa in two communities, a measure of phylogenetic beta diversity
#' @param samp Community data matrix
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param method Character for method used for the swap algorithm ("swap", "tswap", "quasiswap", "backtrack"). If NULL no swap algorithm is applied (uses permatfull from vegan). If mtype="count" the "quasiswap", "swap", "swsh" and "abuswap" methods are available (see details).
#' @param fixedmar Character, stating which of the row/column sums should be preserved ("none", "rows", "columns", "both").
#' @param shuffle Character, indicating whether individuals ("ind"), samples ("samp") or both ("both") should be shuffled, see details.
#' @param strata Numeric vector or factor with length same as nrow(m) for grouping rows within strata for restricted permutations. Unique values or levels are used.
#' @param mtype Matrix data type, either "count" for count data, or "prab" for presence-absence type incidence data.
#' @param burnin Number of null communities discarded before proper analysis in sequential ("swap", "tswap") methods.
#' @param thin Number of discarded permuted matrices between two evaluations in sequential ("swap", "tswap") methods.
#' @param abundance.weighted Should mean nearest taxon distances for each species be weighted by species abundance? (default = FALSE)
#' @param runs Number of randomizations
#' @details See permat (vegan) for detailed options
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
#' @import vegan
#' @export

ses.comdist2 <- function (samp, dis, method = "quasiswap", fixedmar = "both", shuffle = "both", strata = NULL, mtype = "count", burnin = 0, thin = 1, 
                          abundance.weighted = FALSE, runs = 999){
  dis <- as.matrix(dis)
  comdist.obs <- as.matrix(comdist(samp, dis, abundance.weighted = abundance.weighted))

  if(is.null(method)) {
    comdist.rand <- replicate(runs, as.matrix(comdist(permatfull(samp, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, times = 1)$perm[[1]], dis, abundance.weighted)), simplify = FALSE)
  } else {
    comdist.rand <- replicate(runs, as.matrix(comdist(permatswap(samp, method = method, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, burnin = burnin, thin = thin, times = 1)$perm[[1]], dis, abundance.weighted)), simplify = FALSE)
    }
  
  comdist.rand.mean <- apply(X = simplify2array(comdist.rand), MARGIN = 1:2, FUN = mean, na.rm = TRUE)
  
  comdist.rand.sd <- apply(X = simplify2array(comdist.rand), MARGIN = 1:2, FUN = sd, na.rm = TRUE)
  
  comdist.obs.z <- (comdist.obs - comdist.rand.mean)/comdist.rand.sd
  
  comdist.obs.rank <- apply(X = simplify2array(c(list(comdist.obs),comdist.rand)), MARGIN = 1:2, FUN = rank)[1,,]
  comdist.obs.rank <- ifelse(is.na(comdist.rand.mean), NA, comdist.obs.rank)
  diag(comdist.obs.rank) <- NA
  
  comdist.obs.p <- comdist.obs.rank/(runs + 1)
  
  list(ntaxa = specnumber(samp), comdist.obs = comdist.obs, comdist.rand.mean = comdist.rand.mean, 
       comdist.rand.sd = comdist.rand.sd, comdist.obs.rank = comdist.obs.rank, comdist.obs.z = comdist.obs.z, comdist.obs.p = comdist.obs.p, runs = runs)
}
