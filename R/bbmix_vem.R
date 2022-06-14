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
