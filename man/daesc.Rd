% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/daesc.R
\name{daesc}
\alias{daesc}
\title{DAESC for single-cell differential allele-specific expression analysis}
\usage{
daesc(
  y,
  n,
  subj,
  x,
  xnull,
  niter = 200,
  niter_laplace = 2,
  num.nodes = 3,
  optim.method = "BFGS",
  converge_tol = 1e-07
)
}
\arguments{
\item{y}{Alternative allele/haplotype read counts for one SNP/gene. A vector of length equal to the number of cells.}

\item{n}{Total allele-specific read counts. A vector of same length as \code{y}.}

\item{subj}{Donor ID. A vector of same length as \code{y}.}

\item{x}{Design matrix for differential ASE. The number of rows should be equal to \code{length(y)}. Should include an intercept (recommended to be the first column) and variables of interest.}

\item{xnull}{Design matrix under the null hypothesis. Hypothesis testing is based on likelihood ratio test comparing the full model based on \code{x} and null model based on \code{xnull}. The number of rows should be equal to \code{length(y)}.}

\item{niter}{Maximum number of iterations of the variational EM (VEM) algorithm. Default to 200.}

\item{niter_laplace}{Number of Newton-Raphson iterations for estimating the individual-specific random effects at each VEM iteration. This is part of the algorithm for numerical integration. Default to 2. Increasing \code{niter_laplace} leads to more accurate results but slower algorithm.}

\item{num.nodes}{Number of nodes in Gaussian-Hermite quadrature for numeric integration. Default to 3.}

\item{optim.method}{Method of numerical optimization. Can be \code{"BFGS"} or \code{"Nelder-Mead"}.}

\item{converge_tol}{Convergence criterion. VEM algorithm is stopped when the relative increase in log-likelihood to last iteration is less than \code{converge_tol}.}
}
\value{
A list including
\item{b}{Estimate of coefficients representing ASE and differential ASE effects. Each element is log odds ratio of allelic fraction per unit increase in the corresponding column of \code{x}.}
\item{sigma2}{Estimated variance of individual-specific random effects.}
\item{phi}{Estimated over-dispersion parameter in the beta-binomail distribution.}
\item{p}{Mixture probabilities pi0 and 1-pi0}
\item{p.value}{P-value for differential ASE.}
\item{wt}{Posterior probabilities for each individual to be classified into cluster 1 (first column) or cluster 2 (second column).}
\item{llkl}{Log-likelihood.}
\item{llkl.null}{Log-likelihood of the null model.}
\item{note}{Note on convergence status.}
\item{note.null}{Note on convergence status of null model.}
\item{nobs}{Number of cells.}
\item{nsubj}{Number of individuals.}
\item{iter}{Total number of VEM iterations.}
\item{model}{"DAESC-BB" or "DAESC-Mix".}
}
\description{
This function conducts automatic model selection between DAESC-BB (\code{daesc_bb}) and DAESC-Mix (\code{daesc_mix}) based on the number of donors (N). When N<20, \code{daesc} chooses DAESC-BB; when N>=20, \code{daesc} chooses DAESC-Mix.
}
\examples{
library(DAESC)
data("example", package="DAESC")
res <- daesc(y=df$y, n=df$n, subj=df$subj, x=cbind(1,df$x), xnull=matrix(1,nrow=nrow(df),ncol=1), niter=200, niter_laplace=2, num.nodes=3,
optim.method="BFGS", converge_tol=1e-8)
str(res)
}
