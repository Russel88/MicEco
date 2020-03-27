#' Permutation test of Wd* - robust distance-based multivariate analysis of variance
#'
#' Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)
#' 
#' This code is taken from https://github.com/alekseyenko/WdStar/
#' @param dm Distance matrix
#' @param f Factor
#' @param nrep Number of permutations
#' @param strata Factor for permuting in strata
#' @keywords distance
#' @return A list with a p-value, test statitic, and number of permutation
#' @export

WdS.test = function(dm, f, nrep=999, strata=NULL){
    if(!is.factor(f)){
        stop("f has to be a factor")
    }
    if(class(dm) != "dist"){
        stop("dm has to be a distance matrix")
    }
    generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata=strata)
}

WdS = function(dm, f){
    # This method computes Wd* statistic for distance matrix dm and factor f
    ns = table(f)
    dist.ss2 = function(dm2, f){ #dm2 is matrix of square distances; f factor
        K = sapply(levels(f), function(lev) f==lev)
        t(K)%*%dm2%*%K/2
    }
    SS2 = dist.ss2(as.matrix(dm)^2, f)
    s2 = diag(SS2)/ns/(ns-1)
    W = sum(ns/s2)
    
    idxs = apply(combn(levels(f), 2),2, function(idx) levels(f) %in% idx)
    
    Ws = sum(apply(idxs, 2, 
                   function(idx) sum(ns[idx])/prod(s2[idx]) * 
                       (sum(SS2[idx, idx])/sum(ns[idx]) - sum(diag(SS2[idx, idx])/ns[idx]))))
    k=nlevels(f)
    h = sum( (1-ns/s2/W)^2/(ns-1))
    Ws/W/(k-1)/(1+(2*(k-2)/(k^2-1))*h)
}

generic.distance.permutation.test = 
    function(test.statistic, dm, f, nrep=999, strata = NULL){
        N = length(f)
        generate.permutation=function(){
            f[sample(N)]
        }
        
        if(!is.null(strata)){
            # map elements of each strata back to their positions in the factor variable
            strata.map = order(unlist(tapply(seq_along(f), strata, identity)))
            generate.permutation=function(){
                p = unlist(tapply(f,strata,sample)) # permute within strata
                p[strata.map]
            }  
        }
        
        stats = c(test.statistic(dm, f), 
                  replicate(nrep, 
                            test.statistic(dm, generate.permutation())))
        
        p.value = sum(stats>=stats[1])/(nrep+1)
        statistic = stats[1]
        list(p.value = p.value, statistic = statistic, nrep=nrep)
    }

