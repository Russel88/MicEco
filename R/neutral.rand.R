#' Fit Sloan et al. (2006) Neutral Model several times
#'
#' Fit neutral model developed by Sloan et al. (2006, Environ Microbiol 8(4):732-740) and implemented by Burns et al. (2015, ISME J 10(3):655-664) several times on ramdomly picked samples and with 16S rRNA gene copy number corrected rarefaction.
#' @param data A phyloseq object
#' @param n Integer. Number of times to repeat analysis
#' @param s Integer. Number of random samples to for each repetition.
#' @param rRNA What Greengenes database version was used to find OTUs. Atm only "v13.5" is available. Alternatively, A dataframe with two variables: "ID" is the OTU id matched by names in x and "Copy" is the copy number.
#' @param rn Integer. Number of reads to sample for rarefaction
#' @param cores Integer. Number of cores to use for parallel computing.
#' @param naming Optional. A list for naming the output e.g. list(Time="1 Week",Type="Gut")
#' @keywords neutral
#' @return A list of length two; first element contains fit statistics, the second element contains predictions.
#' @import bbmle Hmisc
#' @export

neutral.rand <- function(data,n=NULL,s=NULL,rRNA=NULL,rn=NULL,cores=1,naming=NULL){
  
  # Otu table
  otu <- t(otu_table(data))
  
  # Tax table
  taxa <- tax_table(data)
  
  # Start parallel
  pb <- txtProgressBar(max = n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  results <- foreach(i = 1:n,.packages = c("MicEco","bbmle","Hmisc"),.export = "neutral.fit",.options.snow = opts) %dopar% {
    
    # Subset randomly
    samp <-  otu[sample(nrow(otu),s,replace=FALSE), ]
    
    # Rarefy
    samp.r <- rarefy_rrna(samp,rn,rRNA)
    
    # Remove empty OTUs
    sampx <- samp.r[,colSums(samp.r) > 0]
    
    # Fit model
    fit <- neutral.fit(sampx)
    
    # Insert run number
    stats.sub <- fit[[1]]
    stats.sub$run <- i
    fit[[1]] <- stats.sub
    
    # Insert tax table
    pred.sub <- fit[[2]]
    pred.sub$run <- i
    pred.t <- merge(pred.sub, taxa, by=0, all.x=TRUE)
    pred.t$OTU <- pred.t[,1]
    pred.t <- pred.t[,-1]
    fit[[2]] <- pred.t
    
    return(fit)
  }
  stopCluster(cl)
  close(pb)
  
  # Bind stat dataframe
  results.stats <- as.data.frame(do.call("rbind",lapply(results,function(x) x[[1]])))
  results.stats <- as.data.frame(lapply(results.stats, as.numeric))
  colnames(results.stats) <- c("m","LogLik","gRsqr","N","Samples","OTUs","detection","run")
  
  # Naming
  if(!is.null(naming)){
    results.stats.x <- results.stats
    
    for(i in 1:length(naming)){
      results.stats.x <- cbind(results.stats.x,NA)
      results.stats.x[,ncol(results.stats)+i] <- naming[[i]]
    }
    colnames(results.stats.x) <- c(colnames(results.stats),names(naming))
  } else results.stats.x <- results.stats
  
  # Bind pred dataframe
  results.pred <- as.data.frame(do.call("rbind",lapply(results,function(x) x[[2]])))
  
  # Naming
  if(!is.null(naming)){
    results.pred.x <- results.pred
    
    for(i in 1:length(naming)){
      results.pred.x <- cbind(results.pred.x,NA)
      results.pred.x[,ncol(results.pred)+i] <- naming[[i]]
    }
    colnames(results.pred.x) <- c(colnames(results.pred),names(naming))
  } else results.pred.x <- results.pred
  
  df <- list(results.stats.x,results.pred.x)
  
  return(df)
}
