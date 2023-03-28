#' Obtain initial values for DAESC
daesc_init <- function(y, n, subj, x){
    init.glmm <- glmer(cbind(y,n-y)~x-1+(1|subj), family="binomial")
    sigma2 <- VarCorr(init.glmm)$subj[1]
    sigma2 <- ifelse(sigma2<0.001,0.05,sigma2)

    init.bb <- aod::betabin(cbind(y,n-y)~. -1, random=~1, data=data.frame(y,n,x))
    phi <- 1/(1/init.bb@param[length(init.bb@param)]-1)
    phi <- ifelse(phi<0.001,0.05,phi)
    param.init <- c(summary(init.glmm)$coefficients[,1], sigma2, phi)
    names(param.init)[c(length(param.init)-1,length(param.init))] <- c("sigma2","phi")

    return(param.init)
}

#' DAESC-BB model for single-cell differential allele-specific expression analysis
#' @description DAESC-BB is a beta-binomial random-effects model for differential allele-specific expression (ASE) analysis using multiple single-cell RNA-seq samples. It can perform differential ASE analysis with respect to any conditions represented by the design matrix, regardless of whether the comparison is between cell-types within an individual or across individuals, and regardless of whether the condition of interest is continuous or discrete.
#' @param y Alternative allele/haplotype read counts for one SNP/gene. A vector of length equal to the number of cells.
#' @param n Total allele-specific read counts. A vector of same length as \code{y}.
#' @param subj Donor ID. A vector of same length as \code{y}.
#' @param x Design matrix for differential ASE. The number of rows should be equal to \code{length(y)}.
#' @param xnull Design matrix under the null hypothesis. Hypothesis testing is based on likelihood ratio test comparing the full model based on \code{x} and null model based on \code{xnull}. The number of rows should be equal to \code{length(y)}.
#' @param niter Maximum number of iterations of the variational EM (VEM) algorithm. Default to 200.
#' @param niter_laplace Number of Newton-Raphson iterations for estimating the individual-specific random effects at each VEM iteration. This is part of the algorithm for numerical integration. Default to 2. Increasing \code{niter_laplace} leads to more accurate results but slower algorithm.
#' @param num.nodes Number of nodes in Gaussian-Hermite quadrature for numeric integration. Default to 3.
#' @param optim.method Method of numerical optimization. Can be \code{"BFGS"} or \code{"Nelder-Mead"}.
#' @param converge_tol Convergence criterion. VEM algorithm is stopped when the relative increase in log-likelihood to last iteration is less than \code{converge_tol}.
#'
#' @return A list including
#' \item{b}{Estimate of coefficients representing ASE and differential ASE effects. The first element is the intercept. The following elements are log odds ratios of allelic fraction per unit increase in the variable(s) in \code{x}.}
#' \item{sigma2}{Estimated variance of individual-specific random effects.}
#' \item{phi}{Estimated over-dispersion parameter in the beta-binomail distribution.}
#' \item{p.value}{P-value for differential ASE.}
#' \item{llkl}{Log-likelihood.}
#' \item{llkl.null}{Log-likelihood of the null model.}
#' \item{note}{Note on convergence status.}
#' \item{note.null}{Note on convergence status of null model.}
#' \item{nobs}{Number of cells.}
#' \item{nsubj}{Number of individuals.}
#' \item{iter}{Total number of VEM iterations.}
#'
#' @export
daesc_bb <- function(y, n, subj, x, xnull=NULL, niter=200, niter_laplace=2, num.nodes=3,
                     optim.method="BFGS", converge_tol=1e-7){

    # Convert design matrix to matrix if input is a vector
    if (class(x)[1]!="matrix"){
        X <- matrix(x, ncol=1)
    } else{
        X <- x
    }
    # Initializing
    param.init <- daesc_init(y=y, n=n, subj=subj, x=x)

    print("Fitting DAESC-BB")
    res <- bbmix_vem(param=param.init, y=y, n=n, X=X,
                     subj=subj, niter=niter, niter_laplace=niter_laplace,
                     num.nodes=num.nodes, converge_tol=converge_tol)

    # Convert null design matrix to matrix if input is NULL or a vector
    if (is.null(xnull)){
        XNull <- matrix(1,nrow=length(y),ncol=1)
    } else if (class(x)[1]!="matrix"){
        XNull <- matrix(xnull, ncol=1)
    } else{
        XNull <- xnull
    }
    # Initialize null model
    param.init.null <- daesc_init(y=y, n=n, subj=subj, x=XNull)

    print("Fitting null model")
    res.null <- bbmix_vem(param=param.init.null, y=y, n=n, X=XNull,
                          subj=subj, niter=niter, niter_laplace=niter_laplace,
                          num.nodes=num.nodes, converge_tol=converge_tol)
    res$llkl.null <- res.null$llkl
    res$note.null <- res.null$note
    res$p.value <- pchisq(2*(res$llkl-res$llkl.null),df=ncol(X)-ncol(XNull),lower.tail=F)

    # Clean up results presentation
    names(res$b) <- c("Intercept",paste0("x",1:(length(res$b)-1)))
    res$phi <- as.numeric(res$phi)

    res <- res[c("b","sigma2","phi","p.value","llkl","llkl.null","note","note.null","nobs","nsubj","iter")]
    return(res)
}

#' DAESC-Mix model for single-cell differential allele-specific expression analysis
#' @description DAESC-Mix is a beta-binomial mixture model for differential allele-specific expression (ASE) analysis using multiple single-cell RNA-seq samples. It can perform differential ASE analysis with respect to any conditions represented by the design matrix, regardless of whether the comparison is between cell-types within an individual or across individuals, and regardless of whether the condition of interest is continuous or discrete. DAESC-Mix is an extension of DAESC-BB incorporating implicit haplotype phasing.
#' @param y Alternative allele/haplotype read counts for one SNP/gene. A vector of length equal to the number of cells.
#' @param n Total allele-specific read counts. A vector of same length as \code{y}.
#' @param subj Donor ID. A vector of same length as \code{y}.
#' @param x Design matrix for differential ASE. The number of rows should be equal to \code{length(y)}.
#' @param xnull Design matrix under the null hypothesis. Hypothesis testing is based on likelihood ratio test comparing the full model based on \code{x} and null model based on \code{xnull}. The number of rows should be equal to \code{length(y)}.
#' @param niter Maximum number of iterations of the variational EM (VEM) algorithm. Default to 200.
#' @param niter_laplace Number of Newton-Raphson iterations for estimating the individual-specific random effects at each VEM iteration. This is part of the algorithm for numerical integration. Default to 2. Increasing \code{niter_laplace} leads to more accurate results but slower algorithm.
#' @param num.nodes Number of nodes in Gaussian-Hermite quadrature for numeric integration. Default to 3.
#' @param optim.method Method of numerical optimization. Can be \code{"BFGS"} or \code{"Nelder-Mead"}.
#' @param converge_tol Convergence criterion. VEM algorithm is stopped when the relative increase in log-likelihood to last iteration is less than \code{converge_tol}.
#'
#' @return A list including
#' \item{b}{Estimate of coefficients representing ASE and differential ASE effects. The first element is the intercept. The following elements are log odds ratios of allelic fraction per unit increase in the variable(s) in \code{x}.}
#' \item{sigma2}{Estimated variance of individual-specific random effects.}
#' \item{phi}{Estimated over-dispersion parameter in the beta-binomail distribution.}
#' \item{p}{Mixture probabilities pi0 and 1-pi0}
#' \item{p.value}{P-value for differential ASE.}
#' \item{wt}{Posterior probabilities for each individual to be classified into cluster 1 (first column) or cluster 2 (second column).}
#' \item{llkl}{Log-likelihood.}
#' \item{llkl.null}{Log-likelihood of the null model.}
#' \item{note}{Note on convergence status.}
#' \item{note.null}{Note on convergence status of null model.}
#' \item{nobs}{Number of cells.}
#' \item{nsubj}{Number of individuals.}
#' \item{iter}{Total number of VEM iterations.}
#' @export
daesc_mix <- function(y, n, subj, x, xnull, niter=200, niter_laplace=2, num.nodes=3,
                     optim.method="BFGS", converge_tol=1e-7){
    # param.init <- daesc_init(y=y, n=n, subj=subj, x=x)

    # Convert design matrix to matrix if input is a vector
    if (class(x)[1]!="matrix"){
        X <- matrix(x, ncol=1)
    } else{
        X <- x
    }
    # Initializing
    param.init <- daesc_init(y=y, n=n, subj=subj, x=x)

    print("Fitting DAESC-Mix")
    res <- bbmixture_vem(param=param.init, y=y, n=n, X=X,
                     subj=subj, niter=niter, niter_laplace=niter_laplace,
                     num.nodes=num.nodes, converge_tol=converge_tol)

    # Convert null design matrix to matrix if input is NULL or a vector
    if (is.null(xnull)){
        XNull <- matrix(1,nrow=length(y),ncol=1)
    } else if (class(x)[1]!="matrix"){
        XNull <- matrix(xnull, ncol=1)
    } else{
        XNull <- xnull
    }
    # Initialize null model
    param.init.null <- daesc_init(y=y, n=n, subj=subj, x=XNull)

    print("Fitting null model")
    res.null <- bbmixture_vem(param=param.init.null, y=y, n=n, X=XNull,
                          subj=subj, niter=niter, niter_laplace=niter_laplace,
                          num.nodes=num.nodes, converge_tol=converge_tol)
    res$llkl.null <- res.null$llkl
    res$note.null <- res.null$note
    res$p.value <- pchisq(2*(res$llkl-res$llkl.null),df=ncol(X)-ncol(XNull),lower.tail=F)

    # Clean up results presentation
    names(res$b) <- c("Intercept",paste0("x",1:(length(res$b)-1)))
    res$phi <- as.numeric(res$phi)
    names(res$p) <- c("z=1","z=-1")
    colnames(res$wt) <- c("z=1","z=-1")

    res <- res[c("b","sigma2","phi","p","p.value","wt","llkl","llkl.null","note","note.null","nobs","nsubj","iter")]

    return(res)
}
