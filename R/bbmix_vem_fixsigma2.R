bbmix_vem_fixsigma2 <- function(param, y, n, X, subj, niter, niter_laplace, num.nodes,
                      optim.method="BFGS", converge_tol=1e-6){ # b, sigma2, phi
    b <- param[1:(length(param)-2)]
    sigma2 <- param[length(param)-1]
    phi <- param[length(param)]

    subjfac <- as.integer(factor(subj))
    nsubj <- max(subjfac)

    # Generate Gaussian quadrature
    ghq <- statmod::gauss.quad(num.nodes,kind="hermite")
    ghq$nodes <- ghq$nodes*sqrt(2)
    ghq$weights <- ghq$weights/sqrt(pi)

    # VEM iterations
    llkl <- 10
    randint <- rep(0,nsubj)
    for (iter in 1:niter){
        ## E-step: use normal distribution to approximate the posterior of randint
        linpred <- as.vector(X%*%b)
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
        # The order of randint corresponds to subjfac=1:nsubj (sorted)

        ## M-step
        # Update beta and phi
        bphi <- optim(par=c(b,phi), fn=bbmix_vem_qfunc, control = list(fnscale=-1),
                      method = optim.method, hessian = TRUE,
                      y=y, n=n, X=X, subjfac=subjfac,
                      randint=randint, randint.prec=randint.prec, ghq=ghq)$par
        b <- bphi[1:(length(bphi)-1)]
        phi <- bphi[length(bphi)]

        # Convergence check
        if (iter%%5==0){
            llkl.new <- bbmix_lkl_agq_fixvar(b, sigma2, phi, y, n, X, subj,
                                             niter_laplace=3, num.nodes=1, sumlkl=TRUE)
            print.out <- round(c(iter,b,sigma2,phi,llkl.new),3)
            names(print.out) <- c("iter",paste0("b",1:length(b)),"sigma2","phi","llkl")
            print(print.out)

            converge <- abs(llkl.new-llkl)/abs(llkl)<converge_tol
            if (converge){
                break
            } else{
                llkl <- llkl.new
            }
        }
    }

    if (iter==niter & !converge){
        warning("Algorithm did not converge")
    }

    # # Standard error
    # estimate <- bbmix_vem_se(b,sigma2,phi,y,n,X,subj)
    # return(list(estimate=estimate, llkl=llkl, nobs=length(y), nsubj=nsubj))

    return(list(b=b,sigma2=sigma2,phi=phi,llkl=llkl, nobs=length(y), nsubj=nsubj))
}

## Log-Likelihood after removing the choose(n,y) terms (\sqrt(2*pi) terms are effectively preserved)

