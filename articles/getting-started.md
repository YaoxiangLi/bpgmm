# Getting started with bpgmm

`bpgmm` fits Bayesian parsimonious Gaussian mixture models for
model-based clustering. It targets three posterior inference goals
described by Lu, Li, and Love (2021): the partition of observations, the
number of clusters, and the cluster covariance structure.

## Fit a model

Input data should be a numeric matrix with variables in rows and
observations in columns.

``` r

library(bpgmm)

fit <- pgmm_rjmcmc(
  X = X,
  mInit = 2,
  mVec = c(1, 6),
  qnew = 2,
  burn = 100,
  niter = 1000,
  Mstep = 1,
  Vstep = 1
)
```

Use `Mstep = 1` to allow RJMCMC updates for the number of clusters. Use
`Vstep = 1` to allow RJMCMC updates across the covariance-constraint
models.

## Summarize posterior samples

``` r

summary <- summarize_pgmm_rjmcmc(fit)

summary$Zalloc
summary$nCluster
summary$nConstraint
```

If a reference partition is available, pass it as `trueCluster` to
calculate the adjusted Rand index.

``` r

summary <- summarize_pgmm_rjmcmc(fit, trueCluster = known_labels)
summary$ari
```

## Work with covariance constraints

The paper uses the eight PGMM covariance model labels `CCC`, `CCU`,
`CUC`, `CUU`, `UCC`, `UCU`, `UUC`, and `UUU`. The package stores these
internally as three-number constraint vectors.

``` r

model_to_constraint("UUU")
constraint_to_model(c(1, 0, 0))
```

## Deprecated names

The old names
[`pgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md),
[`summarizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md),
and
[`summerizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md)
remain available as compatibility wrappers, but new code should use
[`pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md)
and
[`summarize_pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md).

## Citation

If you use `bpgmm` in published work, please cite the package and the
methodology paper:

``` r

citation("bpgmm")
```

Lu, X., Li, Y., & Love, T. (2021). On Bayesian Analysis of Parsimonious
Gaussian Mixture Models. *Journal of Classification*, 38, 576-593.
<https://doi.org/10.1007/s00357-021-09391-8>
