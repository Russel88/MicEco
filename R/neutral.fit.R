#' Fit Sloan et al. (2006) Neutral Model
#'
#' Fit neutral model developed by Sloan et al. (2006, Environ Microbiol 8(4):732-740) and implemented by Burns et al. (2015, ISME J 10(3):655-664).
#' @param otu An OTU-table with taxa as columns and samples as rows.
#' @keywords neutral
#' @return A list of length two; first element contains fit statistics, the second element contains predictions.
#' @import bbmle
#' @export

neutral.fit <- function(otu){
  
  # Frequency
  otu.pa <- (otu>0)*1
  freq <- apply(otu.pa, 2, mean)
  
  # Individuals per community
  N <- mean(apply(otu, 1, sum))
  
  # Relative abundance
  p <- apply(otu, 2, function(x) mean(x))/N
  
  # Detection limit
  d = 1/N
  
  # Define likelihood function
  neutral.ll <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    -sum(dnorm(R, 0, sigma,log=TRUE))
  }
  
  # Fit neutral model
  m.mle <- mle2(neutral.ll, start=list(m=0.01, sigma=0.1),method="Nelder-Mead")
  
  # R-squared
  gRsqr <- 1 - exp(-as.numeric(logLik(m.mle))/length(p))
  
  # Predictions
  freq.pred <- pbeta(d, N*m.mle@coef['m']*p, N*m.mle@coef['m']*(1-p), lower.tail=FALSE)
  pred.ci <- binconf(freq.pred*nrow(otu), nrow(otu), return.df=TRUE)
  
  # Bind results
  stats <- cbind(m.mle@coef['m'], m.mle@details$value,gRsqr , N, nrow(otu), length(p), d)
  pred <- cbind(p, freq, freq.pred, pred.ci[,2:3])
  
  results <- list(stats,pred)
  return(results)
}
