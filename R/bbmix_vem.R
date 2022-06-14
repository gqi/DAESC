#' VEM Q function - only terms which are functions of b and phi are preserved
bbmix_vem_qfunc <- function(param, sigma2, y, n, X, subjfac,
                            randint, randint.prec, ghq){
    b <- param[1:(length(param)-1)]
    phi <- param[length(param)]
    linpred <- as.vector(X%*%b)

    zmat <- matrix(rep(ghq$nodes,each=length(linpred)), ncol=length(ghq$nodes))
    mu <- 1/(1+exp(-(linpred+randint[subjfac]+zmat/sqrt(randint.prec)[subjfac])))
    mu[mu<1e-6] <- 1e-6
    mu[mu>1-1e-6] <- 1-1e-6

    # phi = 1/(alpha1+alpha_2); mu = alpha1/(alpha1+alpha2)
    # Warnings are generated here when phi<0
    alpha1 <- mu/phi
    alpha2 <- (1-mu)/phi

    return(sum(colSums(lbeta(alpha1+y,alpha2+n-y)-lbeta(alpha1,alpha2))*ghq$weights))
}

#' Beta-binomial model VEM E-step
#' @description  randint, randint.prec are initial values
bbmix_vem_estep <- function(linpred, sigma2, phi, y, n, subjfac,
                            randint, randint.prec, ghq, niter_laplace){
    # Update the latent random effect variable
    # The order of randint corresponds to subjfac=1:nsubj (sorted)
    for (i in 1:niter_laplace){
        df <- compute_randint_deriv(linpred=linpred, randint=randint,
                                    subjfac=subjfac, y=y, n=n, phi=phi)
        # This is the local optimal random effects
        randint <- randint - 0.9*(df$a1-randint/sigma2)/(df$a2-1/sigma2)
    }

    # Adaptive gaussian quadrature
    df <- compute_randint_deriv(linpred=linpred, randint=randint,
                                subjfac=subjfac, y=y, n=n, phi=phi)
    randint.prec <- abs(df$a2-1/sigma2) # Precision of normal approximate distribution

    return(list(randint=randint, randint.prec=randint.prec))
}

#' Beta-binomial model VEM M-step
#' @param param: b and phi
bbmix_vem_mstep <- function(param, y, n, X, subjfac,
                            randint, randint.prec, ghq, optim.method){
    # Update sigma2
    sigma2 <- mean(randint^2+1/randint.prec)
    # Update beta and phi. par=c(b,phi)
    bphi <- optim(par=param, fn=bbmix_vem_qfunc, control = list(fnscale=-1),
                  method = optim.method, hessian = TRUE,
                  y=y, n=n, X=X, subjfac=subjfac,
                  randint=randint, randint.prec=randint.prec, ghq=ghq)$par
    b <- bphi[1:(length(bphi)-1)]
    phi <- bphi[length(bphi)]

    llkl <- bbmix_lkl_agq_fixvar(b=b, sigma2=sigma2, phi=phi,
                                 y=y, n=n, X=X, subjfac=subjfac,
                                 niter_laplace=3, num.nodes=1, sumlkl=TRUE)
    return(list(b=b, sigma2=sigma2, phi=phi, llkl=llkl))
}


#' Fit beta-binomial random-effects model using variational EM algorithm
#' @export
bbmix_vem <- function(param, y, n, X, subj, niter, niter_laplace, num.nodes,
                           optim.method="BFGS", converge_tol=1e-6){ # b, sigma2, phi

    subjfac <- as.integer(factor(subj))
    nsubj <- max(subjfac)

    # Generate Gaussian quadrature
    ghq <- statmod::gauss.quad(num.nodes,kind="hermite")
    ghq$nodes <- ghq$nodes*sqrt(2)
    ghq$weights <- ghq$weights/sqrt(pi)

    # VEM iterations
    randint <- rep(0,nsubj)
    note <- "Did not converge"
    mstep.out <- list(b=param[1:(length(param)-2)], sigma2=param[length(param)-1],
                      phi=param[length(param)], llkl=-1e300)
    estep.out <- list(randint=rep(0,nsubj), randint.prec=rep(1/mstep.out$sigma2,nsubj))

    for (iter in 1:niter){
        linpred <- as.vector(X%*%mstep.out$b)
        # E-step: use normal distribution to approximate the posterior of randint
        estep.out <- bbmix_vem_estep(linpred=linpred, sigma2=mstep.out$sigma2, phi=mstep.out$phi,
                                    y=y, n=n, subjfac=subjfac, randint=estep.out$randint,
                                    randint.prec=estep.out$randint.prec,
                                    ghq=ghq, niter_laplace=niter_laplace)
        # M-step
        mstep.out.new <- bbmix_vem_mstep(param=c(mstep.out$b,mstep.out$phi), y=y, n=n, X=X, subjfac=subjfac,
                                         randint=estep.out$randint, randint.prec=estep.out$randint.prec,
                                         ghq=ghq, optim.method=optim.method)
        # Convergence check
        llkl_reldiff <- (mstep.out.new$llkl-mstep.out$llkl)/abs(mstep.out$llkl)
        if (llkl_reldiff> -1e-7){
            # If the update improves the likelihood, accept. Otherwise reject and exit
            mstep.out <- mstep.out.new
            if (llkl_reldiff<converge_tol){
                note <- "Converged"
                break
            }
        } else{
            note <- "Update rejected"
            break
        }

        if (iter%%10==0){
            print.out <- round(c(iter,mstep.out$b,mstep.out$sigma2,mstep.out$phi,mstep.out$llkl),3)
            names(print.out) <- c("iter",paste0("b",1:length(mstep.out$b)),"sigma2","phi","llkl")
            print(print.out)
        }
    }

    return(c(mstep.out, list(nobs=length(y), nsubj=nsubj, note=note, iter=iter)))
}

#' Mixture model VEM mixture weights
logBB_bysubj <- function(linpred, sigma2, phi, y, n, subjfac,
                            randint, randint.prec, ghq){

    zmat <- matrix(rep(ghq$nodes,each=length(linpred)), ncol=length(ghq$nodes))
    mu <- 1/(1+exp(-(linpred+randint[subjfac]+zmat/sqrt(randint.prec)[subjfac])))
    mu[mu<1e-6] <- 1e-6
    mu[mu>1-1e-6] <- 1-1e-6

    # phi = 1/(alpha1+alpha_2); mu = alpha1/(alpha1+alpha2)
    # Warnings are generated here when phi<0
    alpha1 <- mu/phi
    alpha2 <- (1-mu)/phi

    logBB_integral <- data.table(v1=as.vector((lbeta(alpha1+y,alpha2+n-y)-lbeta(alpha1,alpha2))%*%ghq$weights),
                                 subjfac=subjfac)
    logBB_integral <- logBB_integral[,list(a1=sum(v1)),by=subjfac] %>%
        arrange(subjfac)

    return(logBB_integral$a1)
}

# Mixture model VEM Q function - only terms that are functions of bet and phi are included
bbmixture_vem_qfunc <- function(param, wt, y, n, X, subjfac,
                            randint, randint.prec, ghq){
    b <- param[1:(length(param)-1)]
    phi <- param[length(param)]
    linpred <- as.vector(X%*%b)
    zmat <- matrix(rep(ghq$nodes,each=length(linpred)), ncol=length(ghq$nodes))

    # Phase 1
    mu <- 1/(1+exp(-(linpred+randint[subjfac]+zmat/sqrt(randint.prec)[subjfac])))
    mu[mu<1e-6] <- 1e-6
    mu[mu>1-1e-6] <- 1-1e-6
    alpha1 <- mu/phi # phi = 1/(alpha1+alpha_2); mu = alpha1/(alpha1+alpha2)
    alpha2 <- (1-mu)/phi

    qfunc <- sum(colSums((lbeta(alpha1+y,alpha2+n-y)-lbeta(alpha1,alpha2))*wt[subjfac,1])*ghq$weights)

    # Phase 2
    mu <- 1/(1+exp(-(-linpred+randint[subjfac]+zmat/sqrt(randint.prec)[subjfac])))
    mu[mu<1e-6] <- 1e-6
    mu[mu>1-1e-6] <- 1-1e-6
    alpha1 <- mu/phi # phi = 1/(alpha1+alpha_2); mu = alpha1/(alpha1+alpha2)
    alpha2 <- (1-mu)/phi

    qfunc <- qfunc + sum(colSums((lbeta(alpha1+y,alpha2+n-y)-lbeta(alpha1,alpha2))*wt[subjfac,2])*ghq$weights)

    return(qfunc)
}

#' Mixture model VEM E-step
#' @description  randint, randint.prec are initial values
bbmixture_vem_estep <- function(linpred, sigma2, phi, p, y, n, subjfac,
                                randint, randint.prec, ghq, niter_laplace){

    # Update the distribution of z_i
    # i-th row of wt corresponds to subjfac==i
    for (i in 1:niter_laplace){
        wt <- cbind(log(p[1]) + logBB_bysubj(linpred, sigma2, phi, y, n, subjfac, randint, randint.prec, ghq),
                    log(p[2]) + logBB_bysubj(-linpred, sigma2, phi, y, n, subjfac, randint, randint.prec, ghq))
        wt <- exp(wt-apply(wt,1,max))
        wt <- wt/rowSums(wt)

        # Update the distribution of a
        df1 <- compute_randint_deriv(linpred, randint, subjfac, y, n, phi, include_deriv0=FALSE)
        df2 <- compute_randint_deriv(-linpred, randint, subjfac, y, n, phi, include_deriv0=FALSE)
        randint <- randint - 0.9*(wt[,1]*df1$a1+wt[,2]*df2$a1-randint/sigma2)/
            (wt[,1]*df1$a2+wt[,2]*df2$a2-1/sigma2)
        df1 <- compute_randint_deriv(linpred, randint, subjfac, y, n, phi, include_deriv0=FALSE)
        df2 <- compute_randint_deriv(-linpred, randint, subjfac, y, n, phi, include_deriv0=FALSE)
        randint.prec <- abs(wt[,1]*df1$a2+wt[,2]*df2$a2-1/sigma2) # Precision of normal approximate distribution
    }

    return(list(wt=wt, randint=randint, randint.prec=randint.prec))
}

#' Mixture model VEM M-step
bbmixture_vem_mstep <- function(param, wt, y, n, X, subjfac,
                                randint, randint.prec, ghq, optim.method){
    # Update p
    p <- colMeans(wt)
    p[p<0.001] <- 0.001
    p <- p/sum(p)
    # Update sigma2
    sigma2 <- mean(randint^2+1/randint.prec)
    # Update beta and phi
    bphi <- optim(par=param, fn=bbmixture_vem_qfunc, control = list(fnscale=-1),
                  method = optim.method, hessian = TRUE,
                  wt=wt, y=y, n=n, X=X, subjfac=subjfac,
                  randint=randint, randint.prec=randint.prec, ghq=ghq)$par
    b <- bphi[1:(length(bphi)-1)]
    phi <- bphi[length(bphi)]

    llkl <- bbmixture_lkl_agq(b=b, sigma2=sigma2, phi=phi, p=p,
                              y=y, n=n, X=X, subjfac=subjfac, niter_laplace=3, num.nodes=1)

    return(list(b=b, sigma2=sigma2, phi=phi, p=p))
}

#' Likelihood function of mixture model
#' @export
bbmixture_lkl_agq <- function(b, sigma2, phi, p, y, n, X, subjfac, niter_laplace,
                              num.nodes){
    llkl1 <- bbmix_lkl_agq_fixvar(b=b, sigma2=sigma2, phi=phi, y=y, n=n, X=X,
                                  subjfac=subjfac, niter_laplace,num.nodes, sumlkl=FALSE)
    llkl2 <- bbmix_lkl_agq_fixvar(b=-b, sigma2=sigma2, phi=phi, y=y, n=n, X=X,
                                  subjfac=subjfac, niter_laplace, num.nodes, sumlkl=FALSE)
    llkl.max <- pmax(llkl1,llkl2)

    return(sum(log(p[1]*exp(llkl1-llkl.max)+p[2]*exp(llkl2-llkl.max))+llkl.max))
}

#' Mixture model VEM wrapper function
#' @description param is initial parameters
#' @export
bbmixture_vem <- function(param, y, n, X, subj, niter, niter_laplace, num.nodes,
                      optim.method="BFGS", converge_tol=1e-6){ # b, sigma2, phi

    subjfac <- as.integer(factor(subj))
    nsubj <- max(subjfac)

    # Generate Gaussian quadrature
    ghq <- statmod::gauss.quad(num.nodes,kind="hermite")
    ghq$nodes <- ghq$nodes*sqrt(2)
    ghq$weights <- ghq$weights/sqrt(pi)

    # VEM iterations
    note <- "Did not converge"
    mstep.out <- list(b=param[1:(length(param)-2)], sigma2=param[length(param)-1],
                      phi=param[length(param)], p=c(0.9,0.1), llkl=-1e300)
    estep.out <- list(randint=rep(0,nsubj), randint.prec=rep(1/mstep.out$sigma2,nsubj))
    for (iter in 1:niter){
        # E-step: use normal distribution to approximate the posterior of randint
        linpred <- as.vector(X%*%mstep.out$b)
        estep.out <- bbmixture_vem_estep(linpred, sigma2=mstep.out$sigma2, phi=mstep.out$phi,
                                         p=mstep.out$p, y, n, subjfac,
                                         randint=estep.out$randint, randint.prec=estep.out$randint.prec,
                                         ghq, niter_laplace)
        # M step
        mstep.out.new <- bbmixture_vem_mstep(param=c(mstep.out$b,mstep.out$phi), wt=estep.out$wt, y, n, X, subjfac,
                                        randint=estep.out$randint, randint.prec=estep.out$randint.prec,
                                        ghq, optim.method=optim.method)

        # Convergence check
        llkl_reldiff <- (mstep.out.new$llkl-mstep.out$llkl)/abs(mstep.out$llkl)
        if (llkl_reldiff> -1e-7){
            # If the update improves the likelihood, accept. Otherwise reject and exit
            mstep.out <- mstep.out.new
            if (llkl_reldiff<converge_tol){
                note <- "Converged"
                break
            }
        } else{
            note <- "Update rejected"
            break
        }

        if (iter%%10==0){
            print.out <- round(c(iter,mstep.out$p[1],mstep.out$b,mstep.out$sigma2,mstep.out$phi,mstep.out$llkl),3)
            names(print.out) <- c("iter","p",paste0("b",1:length(mstep.out$b)),"sigma2","phi","llkl")
            print(print.out)
        }
    }

    # Prepare mixture weight matrix
    wt <- estep.out$wt
    subjdf <- data.frame(subj=subj,subjfac=subjfac) %>%
        distinct() %>% arrange(subjfac)
    rownames(wt) <- subjdf$subj

    return(c(mstep.out, list(llkl=llkl, nobs=length(y), nsubj=nsubj, wt=wt, note=note)))
}
