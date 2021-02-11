#' Clean tax_table such that NAs are replaced with names of the most specific known taxonomy prefixed with the rank.
#'
#' Example: An OTU is annotated as Proteobacteria in the Phylum rank, but NAs in the more specific ranks.
#' All annotations below Phylum for this OTU will be replaced with Phylum_Proteobacteria.
#' @param data \code{phyloseq} or \code{tax_table} object.
#' @return Similar to input, but with NAs replaced in the tax_table
#' @export

ps_tax_clean <- function(data){
    
    # Extract
    if(is(data, "phyloseq")){
        tax <- tax_table(data)
    } else if(is(data, "taxonomyTable")) {
        tax <- data
    } else {
        stop("Input should be a phyloseq or tax_table object")
    }
        
    rankn <- rank_names(data)
    
    # Replace NAs
    tax <- t(apply(tax, 1, function(x) 
        if(sum(is.na(x))>0){
            c(x[!is.na(x)], 
              paste(rankn[max(which(!is.na(x)))], rep(x[max(which(!is.na(x)))], sum(is.na(x))), sep="_"))} 
        else{x}))
    
    # Return
    colnames(tax) <- rankn
    
    if(is(data, "phyloseq")){
        tax_table(data) <- tax_table(tax)
        return(data)
    } else if(is(data, "taxonomyTable")) {
        return(tax_table(tax))
    }
    
}
