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

### Example

```
library(DAESC)
data("example", package="DAESC")
# DAESC-BB analysis
res.bb <- daesc_bb(y=df$y, n=df$n, subj=df$subj, x=df$x, niter=200, niter_laplace=2, num.nodes=3,
                optim.method="BFGS", converge_tol=1e-8)
# DAESC-Mix analysis
res.mix <- daesc_mix(y=df$y, n=df$n, subj=df$subj, x=df$x, niter=200, niter_laplace=2, num.nodes=3,
                optim.method="BFGS", converge_tol=1e-8)
```

Type `?daesc_bb` or `?daesc_mix` in R for detailed documentation.

### Reference

If you use this package, please cite

Qi G, Strober BJ, Popp JM, Ji H,  Battle A. Single-cell allele-specific expression analysis reveals dynamic and cell-type-specific regulatory effects. bioRxiv 2022.10.06.511215 (2022) doi:10.1101/2022.10.06.511215.

### Maintainer

Guanghao Qi (gqi1@jhu.edu)
