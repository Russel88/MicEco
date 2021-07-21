#' Make Venn diagram of shared taxa (ASVs, OTUs) across sample groups
#'
#' Make Venn diagram of shared taxa (ASVs, OTUs) across sample groups from a phyloseq object. Overlap can be weighted by relative abundance
#' #' 
#' Any further arguments to this function are passed to the plot.venn function from the eulerr package. This can be used to change colors, fonts, and other graphical parameters.
#' For example: ps_venn(phy, "Time", quantities = list(type=c("percent","counts"), font = 2), labels = list(cex = 2), col = "red", fill = c("red","blue","green"))
#' @param ps A phyloseq object
#' @param group The grouping factor. Should match variable in sample_data(ps)
#' @param fraction The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count.
#' @param weight If TRUE, the overlaps are weighted by abundance
#' @param relative Should abundances be made relative
#' @param plot If TRUE return a plot, if FALSE return a list with shared and unique taxa
#' @param ... Additional arguments
#' @keywords venn diagram
#' @return An venn plot
#' @import phyloseq
#' @import eulerr
#' @importFrom stats aggregate as.formula
#' @export

ps_venn <- function(ps, group, fraction = 0, weight = FALSE, relative = TRUE, plot = TRUE, ...){
    
    if(relative){
        ps <- transform_sample_counts(ps, function(x) x/sum(x))
    }
    
    if(taxa_are_rows(ps)){
        ps_melted <- reshape2::melt(otu_table(ps))
    } else {
        ps_melted <- reshape2::melt(t(otu_table(ps)))
    }
    ps_melted <- merge(ps_melted, sample_data(ps), by.x = "Var2", by.y = "row.names")

    ps_agg <- aggregate(as.formula(paste("value ~ Var1 +",group)), data = ps_melted, function(x) (sum(x > 0)/length(x) >= fraction) * mean(x))
    
    ps_mat <- reshape2::dcast(as.formula(paste("Var1 ~ ",group)), data = ps_agg, value.var = "value")
    
    rownames(ps_mat) <- ps_mat[, 1]
    ps_mat <- ps_mat[, -1]
    ps_mat_bin <- (ps_mat>0)*1
    
    if(plot){
        if(weight){
            df <- eulerr::venn(ps_mat_bin, weights = rowMeans(ps_mat))
        } else {
            df <- eulerr::venn(ps_mat_bin)
        }
        plot(df, ...)
    } else {
        # Find taxa in all combinations
        combis <- lapply(2:ncol(ps_mat), function(k) lapply(lapply(1:(ncol(combn(1:ncol(ps_mat_bin), m = k))),
                                                                   function(y) ps_mat_bin[, combn(1:ncol(ps_mat_bin), m = k)[, y]]),
                                                            function(x) rownames(x[rowSums(x) >= k, , drop=FALSE])))
        
        # Find taxa in singles
        singles <- apply(ps_mat_bin, 2, function(x) names(x[x > 0]))
        
        # Keep only those NOT in the same combination space
        singles <- lapply(seq_along(singles), function(x) setdiff(singles[[x]], do.call(c, singles[-x])))
        combis <- lapply(combis, function(cc) lapply(seq_along(cc), function(x) setdiff(cc[[x]], do.call(c, cc[-x]))))
        
        # Names
        names(singles) <- colnames(ps_mat_bin)
        for(i in 2:ncol(ps_mat)){
            names(combis[[i-1]]) <- apply(combn(colnames(ps_mat_bin), m = i), 2, function(x) paste(x, collapse = "__"))
        }
        
        # Recursively go through combination space from complex to simple to keep only those in unique combinations
        combis <- rev(combis)
        combis_new <- list()
        for(i in seq_along(combis)){
            if(i == 1) {
                combis_new[[i]] <- combis[[i]]
            } else {
                combis_new[[i]] <- lapply(combis[[i]], function(x) setdiff(x, unlist(combis_new)))
            }
        }
        combis_new <- c(singles, unlist(combis_new, recursive = FALSE))
        return(combis_new[sapply(combis_new, function(x) length(x)>0)])
    }
}
