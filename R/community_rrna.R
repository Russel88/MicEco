#' Calculate the community/sample-wise mean 16S rRNA copy numbers
#'
#' Calculate the average 16S rRNA copy number for each sample in an OTU table 
#' @param x A phyloseq object OR an OTU-table with taxa as rows and OTU names as rownames. OTUs should be picked against the Greengenes v13.5 database, unless a another copy number database is provided.
#' @param copy.database What Greengenes database version was used to find OTUs. Atm only "v13.5" is available. Alternatively, A dataframe with two variables: "ID" is the OTU id matched by rownames in x and "Copy" is the copy number.
#' @param weighted Logical. Should the average copy number be weighted by the relative abundances of the OTUs? If not, it is adviced to rarefy the otu-table first if there are large differences in sample reads.
#' @keywords rrna
#' @return A dataframe with an average copy number for each sample
#' @import phyloseq
#' @export

community_rrna <- function(x, copy.database="v13.5", weighted=TRUE){
  
  # Load Copy database
  if(copy.database == "v13.5"){
    rRNA <- read.table("https://raw.githubusercontent.com/Russel88/MicEco/master/data/gg_13_5_16S.tab")
    colnames(rRNA) <- c("ID","Copy")
  } else {
    rRNA <- copy.database 
  }
  
  # OTU table
  if(class(x) == "phyloseq"){
    if(taxa_are_rows(x)) {
      otu <- as.data.frame(otu_table(x))
    } else{
      otu <- as.data.frame(t(otu_table(x)))
    }
  } else {
    otu <- as.data.frame(x)
  }
    
  # Subset copy number table
  CopyTabSub <- rRNA[rRNA$ID %in% rownames(otu),]
  
  # Normalize or binarize
  if(weighted){
    otu.x <- apply(otu,2,function(x) x/sum(x))
  } else {
    otu.x <- 1*(otu>0)
  }
  
  # Merge
  otu.o <- merge(otu.x,CopyTabSub,by.x="row.names",by.y="ID")
  otu.o <- otu.o[,-1]
  
  # Calculate mean 
  Copy <- apply(otu.o[,1:(ncol(otu.o)-1)], 2, function(x) weighted.mean(otu.o[,ncol(otu.o)], x))
  final <- data.frame(Sample=names(Copy),Copy.number=Copy)
  rownames(final) <- NULL
  return(final)
}
