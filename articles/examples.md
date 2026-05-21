# Worked examples

This article shows a complete, runnable `bpgmm` workflow. The example is
small enough to build quickly on CRAN and pkgdown, but it uses the same
functions as a larger analysis.

## Simulate a small clustering problem

[`pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md)
expects a numeric matrix with variables in rows and observations in
columns. The code below creates two compact clusters in two dimensions.

``` r

library(bpgmm)
#> bpgmm 1.1.6 loaded. If you use bpgmm in published work, please cite it with citation("bpgmm").

set.seed(2026)

X <- cbind(
  matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
  matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
)
known_labels <- rep(1:2, each = 4)

dim(X)
#> [1] 2 8
known_labels
#> [1] 1 1 1 1 2 2 2 2
```

## Choose a covariance model

The paper describes eight covariance structures using three-letter PGMM
labels: `CCC`, `CCU`, `CUC`, `CUU`, `UCC`, `UCU`, `UUC`, and `UUU`. Use
[`model_to_constraint()`](https://yaoxiangli.github.io/bpgmm/reference/model_to_constraint.md)
when calling the sampler, and
[`constraint_to_model()`](https://yaoxiangli.github.io/bpgmm/reference/constraint_to_model.md)
when reading stored results.

``` r

constraint <- model_to_constraint("UUU")
constraint
#> [1] 0 0 0

constraint_to_model(constraint)
#> [1] "UUU"
```

## Fit a short RJMCMC chain

This vignette uses only three posterior iterations so it runs quickly.
For research use, increase `burn` and `niter`, and set `Mstep = 1` and
`Vstep = 1` when you want RJMCMC model selection for the number of
clusters and covariance structure.

``` r

fit_log <- capture.output({
  fit <- pgmm_rjmcmc(
    X = X,
    mInit = 2,
    mVec = c(1, 3),
    qnew = 1,
    burn = 1,
    niter = 3,
    constraint = constraint,
    Mstep = 0,
    Vstep = 0,
    verbose = FALSE
  )
})
tail(fit_log, 1)
#> character(0)
```

The returned object stores posterior samples for allocations, covariance
constraints, component parameters, and hyperparameters.

``` r

names(fit)
#>  [1] "taoList"        "psyList"        "MList"          "lambdaList"    
#>  [5] "YList"          "ZmatList"       "constraintList" "alpha1Vec"     
#>  [9] "alpha2Vec"      "bbetaVec"       "clusIndList"
length(fit$ZmatList)
#> [1] 3
```

## Summarize posterior samples

[`summarize_pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md)
returns the modal allocation, posterior counts for the number of
clusters, posterior counts for covariance models, and optionally an
adjusted Rand index.

``` r

fit_summary <- summarize_pgmm_rjmcmc(fit, trueCluster = known_labels)

fit_summary$Zalloc
#> [1] 2 2 2 2 1 1 1 1
fit_summary$nCluster
#> 
#> 2 
#> 3
fit_summary$nConstraint
#> 
#> UUU 
#>   3
fit_summary$ari
#> [1] 1
```

## Scale the example to real data

For a real analysis, keep the same workflow but use larger sampler
settings. The exact values depend on data size and convergence
diagnostics.

``` r

fit <- pgmm_rjmcmc(
  X = your_data_matrix,
  mInit = 2,
  mVec = c(1, 10),
  qnew = 4,
  burn = 5000,
  niter = 15000,
  constraint = model_to_constraint("UUU"),
  Mstep = 1,
  Vstep = 1,
  SCind = 1,
  verbose = FALSE
)

summarize_pgmm_rjmcmc(fit, trueCluster = known_labels)
```
