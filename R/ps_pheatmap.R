#' Make pheatmap directly from a phyloseq object
#'
#' Make pheatmap directly from a phyloseq object, including agglomoration and filtering. 
#' Default color map is viridis with absent taxa as black.
#' @param ps A phyloseq object
#' @param annot_samp Sample variables to annotate with
#' @param annot_taxa Taxa variables to annotate with
#' @param relative If TRUE, abundances total sum scaled. Default TRUE
#' @param log10 If TRUE, log10 transform abundances. Default TRUE
#' @param tax_agg Taxa rank to agglomorate. Default NULL
#' @param order_taxa If TRUE, taxa are ordered from most to least abundant. Default TRUE
#' @param min_samples Minimum number of samples the features should be present in. Default 1
#' @param min_reads Minimum number of total reads the features should have. Default 1
#' @param min_abundance Minimum mean relative abundance features should have. Default 0
#' @param label_rank Taxa rank to label the taxa. If NULL will label by taxa_names(ps). Default NULL
#' @param color Color palette. Default viridis with absent 
#' @param ... Additional arguments to pheatmap function
#' @return A pheatmap
#' @import phyloseq
#' @import pheatmap
#' @export
ps_pheatmap <- function(ps, annot_samp=NULL, annot_taxa=NULL, relative=TRUE, log10=TRUE, tax_agg=NULL, order_taxa=TRUE, 
                        min_samples=1, min_reads=1, min_abundance=0,
                        label_rank=NULL, color=c("black", viridis::viridis(10)), ...){
    
    # Check
    if(any(!annot_samp %in% sample_variables(ps))){
        stop("Sample annotation not found in phyloseq object")
    }
    if(any(!annot_taxa %in% rank_names(ps))){
        stop("Taxa annotation not found in phyloseq object")
    }
    
    # Agglom
    if(!is.null(tax_agg)){
        ps <- tax_glom(ps, tax_agg)
        label_rank <- tax_agg
    }
    
    # Filter
    if(min_samples > 1 | min_reads > 1 | min_abundance > 0){
        try(ps <- ps_prune(ps, min.samples = min_samples, min.reads = min_reads, min.abundance = min_abundance),
            silent = TRUE)
    }
    
    # Transform
    if(relative){
        ps <- transform_sample_counts(ps, function(x) x/sum(x))
    }
    
    # Extract
    if(taxa_are_rows(ps)){
        otu <- otu_table(ps)
    } else {
        otu <- t(otu_table(ps))
    }
    
    # Log transform
    if(log10){
        otu[otu == 0]  <- min(apply(otu, 2, function(x) min(x[x!=0])))
        otu <- log10(otu)
    }
    
    # Start argument list
    arg_lst <- list(otu)
    
    # Order
    if(order_taxa){
        ordering <- order(rowSums(otu), decreasing = TRUE)
        otu <- otu[ordering, ]
        arg_lst <- list(otu, cluster_rows = FALSE)
    } else {
        ordering <- 1:nrow(otu)
    }
    
    # Annotate samples
    if(!is.null(annot_samp)){
        samp <- sample_data(ps)
        samp <- data.frame(samp[, annot_samp, drop=FALSE])
        arg_lst <- append(arg_lst, list(annotation_col = samp))
    }

    # Annotate taxa
    if(!is.null(annot_taxa)){
        tax <- tax_table(ps)
        tax <- data.frame(tax[, annot_taxa, drop=FALSE])
        arg_lst <- append(arg_lst, list(annotation_row = tax))
    }
    
    # Taxa labels
    if(!is.null(label_rank) | !is.null(tax_agg)){
        tax_labs <- tax_table(ps)[, label_rank]
        tax_labs[is.na(tax_labs)] <- "Other"
        tax_labs <- tax_labs[ordering]
        arg_lst <- append(arg_lst, list(labels_row = tax_labs))
    }
    
    # Add color
    arg_lst <- append(arg_lst, list(color=color))   
    
    # Add additional args
    arg_lst <- append(arg_lst, list(...))
    
    # Plot
    do.call(pheatmap::pheatmap, arg_lst)

}

