#' Rarefaction curve on phyloseq object
#' Theoretical richness with vegan's rarefy function
#' 
#' @param physeq phyloseq object
#' @param subsamp Vector of number of reads to subsample
#' @param trim Remove richness estimations from subsamples larger than the library size
#' @param add_sample_data Add sample data to data.frame
#' @import vegan
#' @import utils
#' @importFrom stats na.omit
#' @export
rcurve <- function(physeq, subsamp = 10^c(1:5), trim = TRUE, add_sample_data = TRUE){
    
    otu <- otu_table(physeq)
    if(!taxa_are_rows(physeq)){
        otu <- t(otu)
    }
    colS <- colSums(otu)
    
    pb <- txtProgressBar(min = 0, max = length(subsamp), style = 3)
    rars <- list()
    for(i in seq_along(subsamp)){
        setTxtProgressBar(pb, i)
        rars[[i]] <- vegan::rarefy(otu, sample = subsamp[i], MARGIN = 2)
    }
    mat <- do.call(cbind, rars)
    
    if(trim){
        mat_bool <- sapply(subsamp, function(i) sapply(colS, function(j) i <= j))
        mat_new <- mat * mat_bool
        mat_new[mat_new == 0] <- NA
    } else {
        mat_new <- mat
    }
    colnames(mat_new) <- subsamp
    
    # To a data.frame
    df <- as.data.frame.table(mat_new)
    colnames(df) <- c("Sample", "Reads", "Richness")
    
    if(trim){
        df <- na.omit(df)
    }

    # Add sample data
    if(add_sample_data){
        samp <- sample_data(physeq)
        df2 <- merge(df, samp, by = "Sample", by.y = "row.names")
        return(df2)
    } else {
        return(df)
    }
    
}
