#' Relevel the Sample variable in a psmelted phyloseq object 
#' 
#' Relevel the Sample variable in a psmelted phyloseq object, 
#' such that similar samples are plotted together with ggplot barcharts 
#' 
#' @param psmelted A phyloseq object melted into a data.frame with psmelt
#' @param ... Arguments passed to hclust
#' @importFrom reshape2 dcast
#' @export
ps_refactor <- function(psmelted, ...){
    if(!is(psmelted, "data.frame")){
        stop("Input should be a phyloseq object melted into a data.frame with psmelt")
    }
    mat <- reshape2::dcast(OTU ~ Sample, value.var = "Abundance", data = psmelted)   
    hc <- hclust(dist(t(mat[, -1])), ...)
    psmelted$Sample <- factor(psmelted$Sample, levels = hc$labels[hc$order])
    return(psmelted)
}    