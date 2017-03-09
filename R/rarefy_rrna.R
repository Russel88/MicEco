#' Rarefy and normalize based on 16S rRNA copy numbers
#'
#' Rarefy an otu-table with the probability of the inverse 16S rRNA copy numbers: The result is a standardized AND normalized otu-table.
#' @param x An otu-table with OTUs as columns OR a phyloseq object
#' @param reads Number of reads to sample. 
#' @param copy.database What Greengenes database version was used to find OTUs. Atm only "v13.5" is available. Alternatively, A dataframe with two variables: "ID" is the OTU id matched by names in x and "Copy" is the copy number.
#' @param seed Random seed for sampling.
#' @keywords rarefy normalize rrna
#' @return A rarefied otu-table
#' @import phyloseq
#' @export

rarefy_rrna <- function(x, reads, copy.database="v13.5", seed=NULL){
  if(class(x) == "phyloseq") y <- rarefy_rrna.phyloseq(x,reads,copy.database,seed) else y <- rarefy_rrna.matrix(x,reads,copy.database,seed)
  return(y)
}

#' @rdname rarefy_rrna
#' @export

rarefy_rrna.matrix <- function (x, reads, rRNA, seed=NULL){
  
  if(is.null(seed)) {
    rand.seed <- sample(1000,1)
    set.seed(rand.seed)
    message(paste("Remember to set seed!","Now set to",rand.seed)) }

  # Load Copy database
  if(copy.database == "v13.5"){
    rRNA <- read.table("https://raw.githubusercontent.com/Russel88/MicEco/master/data/gg_13_5_16S.tab")
    colnames(rRNA) <- c("ID","Copy")
  } else {
    rRNA <- copy.database 
  }
    
  # Probabilities
  rrna <- as.data.frame(rRNA[rRNA$ID %in% colnames(x),])
  rownames(rrna) <- rrna$ID
  rrna <- as.data.frame(t(rrna[,2,drop=FALSE]))
  rrna <- rrna[,order(colnames(rrna))]
  rrna.rev <- 1/rrna
  
  x <- as.matrix(x)
  these.reads <- rep(reads, length = nrow(x))
  x <- x[,order(colnames(x))]
  nm <- colnames(x)
  for (i in 1:nrow(x)) {
    if (sum(x[i, ]) <= these.reads[i]) 
      next
    row <- sample(rep(nm, times = x[i, ]), these.reads[i],prob=rep(rrna.rev, times = x[i, ]))
    row <- table(row)
    ind <- names(row)
    x[i, ] <- 0
    x[i, ind] <- row
  }
  x
}

#' @rdname rarefy_rrna
#' @export

rarefy_rrna.phyloseq <- function (x, reads, rRNA, seed=NULL){
  
  if(is.null(seed)) {
    rand.seed <- sample(1000,1)
    set.seed(rand.seed)
    message(paste("Remember to set seed!","Now set to",rand.seed)) }
  
  x2 <- as.matrix(otu_table(x))
  if(taxa_are_rows(x)) x2 <- as.matrix(t(x2))
  
  # Load Copy database
  if(copy.database == "v13.5"){
    rRNA <- read.table("https://raw.githubusercontent.com/Russel88/MicEco/master/data/gg_13_5_16S.tab")
    colnames(rRNA) <- c("ID","Copy")
  } else {
    rRNA <- copy.database 
  }
    
  # Probabilities
  rrna <- as.data.frame(rRNA[rRNA$ID %in% colnames(x2),])
  rownames(rrna) <- rrna$ID
  rrna <- as.data.frame(t(rrna[,2,drop=FALSE]))
  rrna <- rrna[,order(colnames(rrna))]
  rrna.rev <- 1/rrna
  
  these.reads <- rep(reads, length = nrow(x2))
  x2 <- x2[,order(colnames(x2))]
  nm <- colnames(x2)
  for (i in 1:nrow(x2)) {
    if (sum(x2[i, ]) <= these.reads[i]) 
      next
    row <- sample(rep(nm, times = x2[i, ]), these.reads[i],prob=rep(rrna.rev, times = x2[i, ]))
    row <- table(row)
    ind <- names(row)
    x2[i, ] <- 0
    x2[i, ind] <- row
  }
  if(taxa_are_rows(x)) otu_table(x) <- otu_table(t(x2),taxa_are_rows = TRUE) else otu_table(x) <- otu_table(x2,taxa_are_rows = FALSE)
  x
}

