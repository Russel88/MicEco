#' Prune taxa from phyloseq object by their prevalence or abundance
#'
#' @param data \code{phyloseq} object.
#' @param min.samples Minimum number of samples the features should be present in. Default 0
#' @param min.reads Minimum number of total reads the features should have. Default 0
#' @param min.abundance Minimum mean relative abundance features should have. Default 0
#' @return Similar to input, but with features not reaching the criteria given grouped as "Others"
#' @export

ps_prune <- function(data, min.samples = 0, min.reads = 0, min.abundance = 0){
    
    # Extract from phyloseq
    if(class(data) == "phyloseq"){
        loadNamespace("phyloseq")
        count_table <- phyloseq::otu_table(data)
        if(!phyloseq::taxa_are_rows(data)) count_table <- t(count_table)
    } else {
        stop("Input should be a phyloseq object")
    }
    
    # Check
    if(any("Others" %in% rownames(count_table))) stop("There is already a feature named 'Others' in data. Maybe data has been pre-processed already")
    
    # Exclude by samples
    exclude.samples <- which(rowSums(count_table > 0) < min.samples)
    
    # Exclude by reads
    exclude.reads <- which(rowSums(count_table) < min.reads)
    
    # Exclude by abundance
    count_rel <- apply(count_table, 2, function(x) x/sum(x))
    exclude.abundance <- which(rowMeans(count_rel) < min.abundance)
    
    # Which to exclude
    exclude <- unique(c(exclude.samples,exclude.reads,exclude.abundance))
    if(length(exclude) == 0) stop("No features to group!")
    if(length(exclude) == nrow(count_table)) stop("All features would be grouped as 'Others'!")
    message(paste(length(exclude),"features grouped as 'Others' in the output"))
    
    # Make new count_table
    count.keep <- count_table[-exclude,]
    count.other <- count_table[exclude,]
    count.other <- colSums(count.other)
    count.new <- rbind(count.keep,count.other)
    rownames(count.new)[nrow(count.new)] <- "Others"
    
    # Output
    # Fix tax_table
    tax <- as.data.frame(unclass(phyloseq::tax_table(data)))
    tax.keep <- tax[-exclude,]
    tax.new <- rbind(tax.keep,NA)
    rownames(tax.new)[nrow(tax.new)] <- "Others"
    if(phyloseq::taxa_are_rows(data)){
        data.new <- phyloseq::phyloseq(phyloseq::otu_table(count.new, taxa_are_rows = TRUE),
                                       phyloseq::sample_data(data),
                                       phyloseq::tax_table(as.matrix(tax.new)))
    } 
    if(!phyloseq::taxa_are_rows(data)){
        data.new <- phyloseq::phyloseq(phyloseq::otu_table(t(count.new), taxa_are_rows = FALSE),
                                       phyloseq::sample_data(data),
                                       phyloseq::tax_table(as.matrix(tax.new)))
    } 
    return(data.new)
    
}
