#' Fit Sloan et al. (2006) Neutral Model
#'
#' Fit neutral model developed by Sloan et al. (2006, Environ Microbiol 8(4):732-740) and implemented by Burns et al. (2015, ISME J 10(3):655-664).
#' #' Outputs three different pseudo R2 values for model fit:
#' - McFadden's R2: 1 - (logLik(model) / logLik(null model))
#' - Cox & Snell's R2: 1 - exp((2/n) * (logLik(null model) - logLik(model)))
#' - Nagelkerke's R2: Cox & Snell's R2 / (1 - exp((2/n) * logLik(null model)))
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
  
  # extract fitted log-likelihood (numeric)
  LL_m <- as.numeric(logLik(m.mle))

  # build null model: predict = mean(freq), estimate sigma0 from residuals
  pred_null <- mean(freq)
  resid0 <- freq - pred_null
  # MLE sigma for normal residuals is sqrt(mean(resid^2))
  sigma0 <- sqrt(mean(resid0^2))
  # compute null log-likelihood (sum of normal log-densities; n = number of species)
  n <- length(freq)
  LL_0 <- sum(dnorm(resid0, mean=0, sd=sigma0, log=TRUE))
  # (equivalently LL_0 = -0.5 * n * (log(2*pi) + log(sigma0^2) + 1))

  # McFadden
  R2_McF <- 1 - (LL_m / LL_0)

  # Cox & Snell
  R2_CS <- 1 - exp((2 / n) * (LL_0 - LL_m))

  # Simple Nagelkerke commonly uses: R2_N = R2_CS / (1 - exp( (2/n) * LL_0 ))
  R2_Nag <- R2_CS / (1 - exp((2 / n) * LL_0))

  # Predictions
  freq.pred <- pbeta(d, N*m.mle@coef['m']*p, N*m.mle@coef['m']*(1-p), lower.tail=FALSE)
  pred.ci <- binconf(freq.pred*nrow(otu), nrow(otu), return.df=TRUE)
    
  # Bind results
  stats <- cbind(m.mle@coef['m'], m.mle@details$value, R2_McF, R2_CS, R2_Nag, N, nrow(otu), length(p), d)
  pred <- cbind(p, freq, freq.pred, pred.ci[,2:3])
  
  results <- list(stats,pred)
  return(results)
}
