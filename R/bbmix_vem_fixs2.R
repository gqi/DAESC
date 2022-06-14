#' Beta-binomial model VEM M-step
#' @param param: b and phi
bbmix_vem_mstep_fixs2 <- function(param, sigma2, y, n, X, subjfac,
                            randint, randint.prec, ghq, optim.method){
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
bbmix_vem_fixs2 <- function(param, y, n, X, subj, niter, niter_laplace, num.nodes,
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
        mstep.out.new <- bbmix_vem_mstep_fixs2(param=c(mstep.out$b,mstep.out$phi), sigma2=mstep.out$sigma2, y=y, n=n, X=X, subjfac=subjfac,
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
