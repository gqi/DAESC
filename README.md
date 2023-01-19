# DAESC: Differential Allelic Expression Analysis using Single-Cell data

DAESC is a software for differential allele-specific expression (ASE) analysis using single-cell RNA-seq data of multiple individuals. It can be applied to any comparison represented by a design matrix, including but not limited to: 1) discrete cell types 2) continuous cell states (e.g. pseudotime) 3) case-control disease status 4) continuous phenotype of the donors (e.g. BMI, blood pressure). DAESC includes two components: DAESC-BB and DAESC-Mix.

* **DAESC-BB** is a beta-binomial regression model with individual-specific random-effects to capture the sample repeat structure arising from having multiple cells per individual. It is the baseline model of DAESC and can be used regardless of sample size.
* **DAESC-Mix** is an extended model from DAESC-BB incorporating implicit haplotype phasing. Implicit phasing is conducted against the the cis-regulatory genetic variant or epigenetic alteration that causes ASE, which is usually **unobserved**. In the eQTL setting, for example, the expression-increasing allele of the eQTL SNP can be on the haplotype of either the alternative or the reference allele of the exonic SNP where ASE is assessed (eSNP). As a result, different individuals may have opposite allelic imbalance actually representing the same regulatory effect. If not addressed, allelic imbalance will cancel each other across individuals, leading to diminished signal. DAESC-Mix uses latent variables to capture both possibilities, leading to a mixture model. DAESC-Mix can be used when the number of individuals is reasonably large (e.g. N>=20).

### Installation

```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("gqi/DAESC")
``` 

Installation typically takes <30s.

### Example

```
library(DAESC)
data("example", package="DAESC")
# DAESC-BB analysis (run time: ~17s)
res.bb <- daesc_bb(y=df$y, n=df$n, subj=df$subj, x=df$x, niter=200, niter_laplace=2, num.nodes=3,
                optim.method="BFGS", converge_tol=1e-8)
# DAESC-Mix analysis (run time: ~46s)
res.mix <- daesc_mix(y=df$y, n=df$n, subj=df$subj, x=df$x, niter=200, niter_laplace=2, num.nodes=3,
                optim.method="BFGS", converge_tol=1e-8)
```

View results

```
str(res.mix)
#List of 13
# $ b        : Named num [1:2] -0.17 0.979
#  ..- attr(*, "names")= chr [1:2] "Intercept" "x1"
# $ sigma2   : num 0.0564
# $ phi      : num 0.672
# $ p        : Named num [1:2] 0.765 0.235
#  ..- attr(*, "names")= chr [1:2] "z=1" "z=-1"
# $ p.value  : num 2.92e-15
# $ wt       : num [1:53, 1:2] 1.00 4.98e-01 1.03e-01 1.75e-06 9.99e-01 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:53] "HPSI0114i-eipl_1" "HPSI0114i-iisa_1" "HPSI0114i-iisa_3" "HPSI0114i-joxm_1" ...
#  .. ..$ : chr [1:2] "z=1" "z=-1"
# $ llkl     : num -27090
# $ llkl.null: num -27122
# $ note     : chr "Converged"
# $ note.null: chr "Converged"
# $ nobs     : int 3792
# $ nsubj    : int 53
# $ iter     : int 27
```

Type `?daesc_bb` or `?daesc_mix` in R for detailed documentation.

### Reference

If you use this package, please cite

Qi G, Strober BJ, Popp JM, Ji H,  Battle A. Single-cell allele-specific expression analysis reveals dynamic and cell-type-specific regulatory effects. bioRxiv 2022.10.06.511215 (2022) doi:10.1101/2022.10.06.511215.

### Maintainer

Guanghao Qi (gqi1@jhu.edu)

### Dependencies

The software can be run on any operating system as long as R (>= 4.1) is installed. Below is a list of required packages. 
```
data.table,
dplyr,
lme4,
aod,
statmod,
numDeriv
```

DAESC has been tested successfully on R/4.1 and R/4.0.
