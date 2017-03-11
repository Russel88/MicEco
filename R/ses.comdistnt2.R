#' Standardized effect size of inter-community MNTD (betaMNTD, betaNTI)
#'
#' Standardized effect size of MNTD (mean nearest taxon distance) separating taxa in two communities, a measure of phylogenetic beta diversity
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
#' @import vegan
#' @export

ses.comdistnt2 <- function (samp, dis, method = "quasiswap", fixedmar = "both", shuffle = "both", strata = NULL, mtype = "count", burnin = 0, thin = 1, 
                          abundance.weighted = FALSE, exclude.conspecifics = FALSE, runs = 999){
  dis <- as.matrix(dis)
  comdistnt.obs <- as.matrix(comdistnt(samp, dis, abundance.weighted = abundance.weighted, exclude.conspecifics = exclude.conspecifics))

  if(is.null(method)) {
    comdistnt.rand <- replicate(runs, as.matrix(comdistnt(permatfull(samp, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, times = 1)$perm[[1]], dis, abundance.weighted, exclude.conspecifics)), simplify = FALSE)
  } else {
    comdistnt.rand <- replicate(runs, as.matrix(comdistnt(permatswap(samp, method = method, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, burnin = burnin, thin = thin, times = 1)$perm[[1]], dis, abundance.weighted, exclude.conspecifics)), simplify = FALSE)
    }
  
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
