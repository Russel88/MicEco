#' Standardized effect size of unifrac
#'
#' Standardized effect size of unifrac
#' @param physeq phyloseq-class, containing at minimum a phylogenetic tree and otu table
#' @param method "taxa.labels" shuffles labels in phylogenetic tree (Ignores fixedmar, shuffle, strata, mtype). If NULL then no swap algorithm is applied (i.e. uses permatfull from vegan). Else the method used for the swap algorithm ("swap", "tswap", "quasiswap", "backtrack"). If mtype="count" the "quasiswap", "swap", "swsh" and "abuswap" methods are available (see details).
#' @param fixedmar Character, stating which of the row/column sums should be preserved ("none", "rows", "columns", "both").
#' @param shuffle Character, indicating whether individuals ("ind"), samples ("samp") or both ("both") should be shuffled.
#' @param strata Numeric vector or factor for grouping samples within strata for restricted permutations. Unique values or levels are used.
#' @param mtype Matrix data type, either "count" for count data, or "prab" for presence-absence type incidence data.
#' @param burnin Number of null communities discarded before proper analysis in sequential ("swap", "tswap") methods.
#' @param thin Number of discarded permuted matrices between two evaluations in sequential ("swap", "tswap") methods.
#' @param weighted Should unifrac be weighted by species abundance? (default = TRUE)
#' @param normalized (Optional). Logical. Should the output be normalized such that values range from 0 to 1 independent of branch length values? Default is TRUE. Note that (unweighted) UniFrac is always normalized by total branch-length, and so this value is ignored when weighted == FALSE.
#' @param runs Number of randomizations
#' @param cores Number of cores to use for UniFrac of null communities. Default is 1
#' @details See permat (vegan) for detailed options on permutation
#' @keywords ses unifrac
#' @return A list of results:
#' \itemize{
#'   \item unifrac.obs - Observed unifrac between communities
#'   \item unifrac.rand.mean - Mean unifrac between null communities
#'   \item unifrac.rand.sd - Standard deviation of unifrac between null communities
#'   \item unifrac.obs.rank - Rank of observed unifrac vs. null unifrac
#'   \item unifrac.obs.z - Standardized effect size of unifrac vs. null unifrac (= (unifrac.obs - unifrac.rand.mean) / unifrac.rand.sd)
#'   \item unifrac.obs.p - P-value (quantile) of observed unifrac vs. null communities (= unifrac.obs.rank / runs + 1)
#' }
#' @import picante
#' @import vegan
#' @import phyloseq
#' @import foreach
#' @import doParallel
#' @export

ses.UniFrac <- function (physeq, method = "taxa.labels", fixedmar = "both", shuffle = "both", strata = NULL, mtype = "count", burnin = 0, thin = 1, 
                         weighted = TRUE, normalized = TRUE, runs = 99, cores = 1){
  
  # Observed
  unifrac.obs <- as.matrix(UniFrac(physeq, weighted = weighted, normalized = normalized))
  
  # Extract otu-table
  if(taxa_are_rows(physeq)) otu <- as.data.frame(t(otu_table(physeq))) else otu <- as.data.frame(otu_table(physeq))
  if(taxa_are_rows(physeq) & fixedmar == "rows") fixedmar <- "columns" 
  if(taxa_are_rows(physeq) & fixedmar == "columns") fixedmar <- "rows" 
  
  # Create nulls
  if(method == "taxa.labels") {
    rands <- replicate(runs, tipShuffle(phy_tree(physeq)), simplify = FALSE)
  } else {
    if(is.null(method)) {
      rands <- permatfull(otu, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, times = runs)$perm
    } else {
      rands <- permatswap(otu, method = method, fixedmar = fixedmar, shuffle = shuffle, strata = strata, mtype = mtype, burnin = burnin, thin = thin, times = runs)$perm
    }
  }
  
  # Temporary data
  temp <- physeq
  
  # Unifrac of random communites
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  unifrac.rand <- foreach(i = 1:runs, .packages = "phyloseq") %dopar% {
    
    if(method == "taxa.labels"){
      phy_tree(temp) <- phy_tree(rands[[i]])
    } else {
      otu_table(temp) <- otu_table(rands[[i]], taxa_are_rows = FALSE)
    }
    
    unifrac.rand.sub <- as.matrix(UniFrac(temp, weighted = weighted, normalized = normalized))
  }  
  stopCluster(cl)
  
  unifrac.rand.mean <- apply(X = simplify2array(unifrac.rand), MARGIN = 1:2, FUN = mean, na.rm = TRUE)
  
  unifrac.rand.sd <- apply(X = simplify2array(unifrac.rand), MARGIN = 1:2, FUN = sd, na.rm = TRUE)
  
  unifrac.obs.z <- (unifrac.obs - unifrac.rand.mean)/unifrac.rand.sd
  
  unifrac.obs.rank <- apply(X = simplify2array(c(list(unifrac.obs),unifrac.rand)), MARGIN = 1:2, FUN = rank)[1,,]
  unifrac.obs.rank <- ifelse(is.na(unifrac.rand.mean), NA, unifrac.obs.rank)
  diag(unifrac.obs.rank) <- NA
  
  unifrac.obs.p <- unifrac.obs.rank/(runs + 1)
  
  list(unifrac.obs = unifrac.obs, unifrac.rand.mean = unifrac.rand.mean, 
       unifrac.rand.sd = unifrac.rand.sd, unifrac.obs.rank = unifrac.obs.rank, unifrac.obs.z = unifrac.obs.z, unifrac.obs.p = unifrac.obs.p)
}