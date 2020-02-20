#' Make Euler diagram of shared taxa (ASVs, OTUs) across sample groups
#'
#' Make Euler diagram of shared taxa (ASVs, OTUs) across sample groups from a phyloseq object. Overlap can be weighted by relative abundance
#' @param ps A phyloseq object
#' @param group The grouping factor. Should match variable in sample_data(ps)
#' @param fraction The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count.
#' @param weight If TRUE, the overlaps are weighted by abundance
#' @param type "percent" or "counts"
#' @param relative Should abundances be made relative
#' @param ... Additional arguments
#' @keywords euler diagram
#' @return An euler plot
#' @import phyloseq
#' @import eulerr
#' @importFrom stats aggregate as.formula
#' @export

ps_euler <- function(ps, group, fraction = 0, weight = FALSE, type = "percent", relative = TRUE, ...){
    
    if(relative){
        ps <- transform_sample_counts(ps, function(x) x/sum(x))
    }
    
    ps_melted <- reshape2::melt(otu_table(ps))
    ps_melted <- merge(ps_melted, sample_data(ps), by.x = "Var2", by.y = "row.names")
    
    ps_agg <- aggregate(as.formula(paste("value ~ Var1 +",group)), data = ps_melted, function(x) (sum(x > 0)/length(x) >= fraction) * mean(x))
    
    ps_mat <- reshape2::dcast(as.formula(paste("Var1 ~ ",group)), data = ps_agg, value.var = "value")
    
    ps_mat <- ps_mat[, -1]
    ps_mat_bin <- (ps_mat>0)*1
    
    if(weight){
        df <- eulerr::euler(ps_mat_bin, weights = rowMeans(ps_mat))
    } else {
        df <- eulerr::euler(ps_mat_bin)
    }
    
    plot(df, quantities = list(type=type), ...)
    
}
