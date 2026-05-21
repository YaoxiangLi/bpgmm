# Getting started with bpgmm

`bpgmm` fits Bayesian parsimonious Gaussian mixture models for
model-based clustering. It targets three posterior inference goals
described by Lu, Li, and Love (2021): the partition of observations, the
number of clusters, and the cluster covariance structure.

For the model equations, prior distributions, PGMM covariance labels,
and RJMCMC move types, see
[`vignette("model-and-sampler", package = "bpgmm")`](https://yaoxiangli.github.io/bpgmm/articles/model-and-sampler.md).

## Fit a model

Input data should be a numeric matrix with variables in rows and
observations in columns. The small example below has two variables and
eight observations. It is intentionally short so the vignette runs
quickly; use larger `burn` and `niter` values for real analysis.

``` r

library(bpgmm)
#> bpgmm 1.1.6 loaded. If you use bpgmm in published work, please cite it with citation("bpgmm").

set.seed(2026)

X <- cbind(
  matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
  matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
)
known_labels <- rep(1:2, each = 4)

fit_log <- capture.output({
  fit <- pgmm_rjmcmc(
    X = X,
    mInit = 2,
    mVec = c(1, 3),
    qnew = 1,
    burn = 1,
    niter = 3,
    constraint = model_to_constraint("UUU"),
    Mstep = 0,
    Vstep = 0,
    verbose = FALSE
  )
})
tail(fit_log, 1)
#> character(0)
```

The arguments control the sampler:

- `mInit` is the initial number of clusters.
- `mVec` is the allowed range for the number of clusters.
- `qnew` is the number of latent factors for a new cluster.
- `constraint` selects the starting covariance model.
- `Mstep = 1` allows RJMCMC updates for the number of clusters.
- `Vstep = 1` allows RJMCMC updates across covariance-constraint models.

A longer run for applied work might look like this:

``` r

fit <- pgmm_rjmcmc(
  X = methylation_matrix,
  mInit = 2,
  mVec = c(1, 6),
  qnew = 2,
  burn = 100,
  niter = 1000,
  Mstep = 1,
  Vstep = 1,
  verbose = FALSE
)
```

## Summarize posterior samples

``` r

summary <- summarize_pgmm_rjmcmc(fit, trueCluster = known_labels)

summary$Zalloc
#> [1] 2 2 2 2 1 1 1 1
summary$nCluster
#> 
#> 2 
#> 3
summary$nConstraint
#> 
#> UUU 
#>   3
summary$ari
#> [1] 1
```

If a reference partition is available, pass it as `trueCluster` to
calculate the adjusted Rand index.

## Work with covariance constraints

The paper uses the eight PGMM covariance model labels `CCC`, `CCU`,
`CUC`, `CUU`, `UCC`, `UCU`, `UUC`, and `UUU`. The package stores these
internally as three-number constraint vectors.

``` r

model_to_constraint("UUU")
#> [1] 0 0 0
constraint_to_model(c(1, 0, 0))
#> [1] "CUU"
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
#> To cite package 'bpgmm' in publications use:
#> 
#>   Lu X, Li Y, Love T (2021). On Bayesian Analysis of
#>   Parsimonious Gaussian Mixture Models. Journal of
#>   Classification, 38, 576-593. doi:10.1007/s00357-021-09391-8
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {On Bayesian Analysis of Parsimonious Gaussian Mixture Models},
#>     author = {Xiang Lu and Yaoxiang Li and Tanzy Love},
#>     journal = {Journal of Classification},
#>     year = {2021},
#>     volume = {38},
#>     pages = {576--593},
#>     doi = {10.1007/s00357-021-09391-8},
#>   }
```

Lu, X., Li, Y., & Love, T. (2021). On Bayesian Analysis of Parsimonious
Gaussian Mixture Models. *Journal of Classification*, 38, 576-593.
<https://doi.org/10.1007/s00357-021-09391-8>
