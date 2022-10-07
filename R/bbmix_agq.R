#' Compute first and second derivatives of
#' sum(lbeta(mu/phi+y,(1-mu)/phi+n-y)-lbeta(mu/phi,(1-mu)/phi)) wrt random intercept
#' length(randint)==nsubj, the k-th element corresponds to subjfac==k
compute_randint_deriv <- function(linpred, randint, subjfac, y, n, phi, include_deriv0=FALSE){
    mu <- 1/(1+exp(-(linpred + randint[subjfac]))) # linpred should be a vector, not matrix
    mu[mu<0.0000001] <- 0.0000001
    mu[mu>0.9999999] <- 0.9999999

    # Digamma function becomes NaN when the argument is 0
    temp1 <- digamma(mu/phi+y)-digamma((1-mu)/phi+n-y)-
        digamma(mu/phi)+digamma((1-mu)/phi)
    temp2 <- trigamma(mu/phi+y)+trigamma((1-mu)/phi+n-y)-
        trigamma(mu/phi)-trigamma((1-mu)/phi)

    if (include_deriv0){
        df <- data.table(v1=temp1*mu*(1-mu)/phi, v2=temp2*mu^2*(1-mu)^2/phi^2+temp1*(1-2*mu)*mu*(1-mu)/phi,
                         v3=lbeta(mu/phi+y,(1-mu)/phi+n-y)-lbeta(mu/phi,(1-mu)/phi), subjfac=subjfac)
        df <- df[,list(a1=sum(v1),a2=sum(v2),a3=sum(v3)),by=subjfac]

    } else{
        df <- data.table(v1=temp1*mu*(1-mu)/phi, v2=temp2*mu^2*(1-mu)^2/phi^2+temp1*(1-2*mu)*mu*(1-mu)/phi,
                         subjfac=subjfac)
        df <- df[,list(a1=sum(v1),a2=sum(v2)),by=subjfac]
    }

    # By default df$subjfac==unique(subjfac), need to resort df by increasing order of subjfac
    return(df%>%arrange(subjfac))
}

#' Log-Likelihood after removing the choose(n,y) terms (1/sqrt(2*pi) terms are effectively preserved)
#' @param If num.nodes=1, likelihood is approximated by Laplace method. num.nodes>1 is deprecated.
#' @import lme4
#' @import numDeriv
#' @import data.table
bbmix_lkl_agq_fixvar <- function(b, sigma2, phi, y, n, X, subjfac, niter_laplace,
                                 num.nodes, sumlkl=TRUE){
    linpred <- as.vector(X%*%b)
    nsubj <- max(subjfac)
    randint <- rep(0,nsubj)

    for (i in 1:niter_laplace){
        df <- compute_randint_deriv(linpred=linpred, randint=randint,
                                    subjfac=subjfac, y=y, n=n, phi=phi)
        # This is the local optimal random effects
        randint <- randint - 0.9*(df$a1-randint/sigma2)/(df$a2-1/sigma2)
    }

    if (num.nodes==1){
        # Laplace approximation
        df <- compute_randint_deriv(linpred=linpred, randint=randint,
                                    subjfac=subjfac, y=y, n=n, phi=phi, include_deriv0 = TRUE)
        randint.prec <- abs(df$a2-1/sigma2) # Precision of normal approximate distribution

        fval <- df$a3 - randint^2/sigma2/2 # Mf(x0)

        if (sumlkl){
            llkl <- sum(fval-0.5*log(randint.prec))-nsubj/2*log(sigma2)
        } else{
            llkl <- fval-0.5*log(randint.prec)-log(sigma2)/2
        }
    } else{
        stop("num.nodes should be a 1")
    }

    return(llkl)
}

#' Likelihood function of mixture model
bbmixture_lkl_agq <- function(b, sigma2, phi, p, y, n, X, subjfac, niter_laplace,
                              num.nodes){
    llkl1 <- bbmix_lkl_agq_fixvar(b=b, sigma2=sigma2, phi=phi, y=y, n=n, X=X,
                                  subjfac=subjfac, niter_laplace,num.nodes, sumlkl=FALSE)
    llkl2 <- bbmix_lkl_agq_fixvar(b=-b, sigma2=sigma2, phi=phi, y=y, n=n, X=X,
                                  subjfac=subjfac, niter_laplace, num.nodes, sumlkl=FALSE)
    llkl.max <- pmax(llkl1,llkl2)

    return(sum(log(p[1]*exp(llkl1-llkl.max)+p[2]*exp(llkl2-llkl.max))+llkl.max))
}
