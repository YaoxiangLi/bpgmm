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
#> bpgmm 1.2.4 loaded. If you use bpgmm in published work, please cite it with citation("bpgmm").

set.seed(2026)

X <- cbind(
  matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
  matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
)
known_labels <- rep(1:2, each = 4)
```

The example is easiest to check visually. Each point is one observation,
and the colors show the reference labels used later for the ARI
calculation.

``` r

plot(
  X[1, ], X[2, ],
  col = c("#0072B2", "#D55E00")[known_labels],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Small two-cluster example",
  asp = 1
)
legend(
  "topleft",
  legend = paste("Reference", sort(unique(known_labels))),
  col = c("#0072B2", "#D55E00"),
  pch = 19,
  bty = "n"
)
```

![Scatter plot of the small two-cluster
example.](getting-started_files/figure-html/unnamed-chunk-3-1.png)

``` r

fit_log <- capture.output({
  fit <- pgmm_rjmcmc(
    X = X,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 1,
    niter = 3,
    constraint = model_to_constraint("UUU"),
    m_step = 0,
    v_step = 0,
    verbose = FALSE
  )
})
tail(fit_log, 1)
#> character(0)
```

The arguments control the sampler:

- `m_init` is the initial number of clusters.
- `m_range` is the allowed range for the number of clusters.
- `q_new` is the number of latent factors for a new cluster.
- `constraint` selects the starting covariance model.
- `m_step = 1` allows RJMCMC updates for the number of clusters.
- `v_step = 1` allows RJMCMC updates across covariance-constraint
  models.

A longer run for applied work might look like this:

``` r

fit <- pgmm_rjmcmc(
  X = methylation_matrix,
  m_init = 2,
  m_range = c(1, 6),
  q_new = 2,
  burn = 100,
  niter = 1000,
  m_step = 1,
  v_step = 1,
  verbose = FALSE
)
```

## Summarize posterior samples

``` r

summary <- summarize_pgmm_rjmcmc(fit, true_cluster = known_labels)

summary$allocation
#> [1] 2 2 2 2 1 1 1 1
summary$n_clusters
#> 
#> 2 
#> 3
summary$n_constraints
#> 
#> UUU 
#>   3
summary$ari
#> [1] 1
```

If a reference partition is available, pass it as `true_cluster` to
calculate the adjusted Rand index.

The summary allocation can be plotted back on the original coordinates.
This is a useful first diagnostic before moving on to longer chains and
convergence checks.

``` r

plot(
  X[1, ], X[2, ],
  col = c("#009E73", "#CC79A7", "#E69F00")[summary$allocation],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Posterior modal allocation",
  asp = 1
)
text(X[1, ], X[2, ], labels = seq_along(summary$allocation), pos = 3, cex = 0.75)
legend(
  "topleft",
  legend = paste("Cluster", sort(unique(summary$allocation))),
  col = c("#009E73", "#CC79A7", "#E69F00")[sort(unique(summary$allocation))],
  pch = 19,
  bty = "n"
)
```

![Scatter plot colored by posterior modal
allocation.](getting-started_files/figure-html/unnamed-chunk-7-1.png)

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

## Naming convention

Starting with version 1.2.0, the public API uses snake_case names
throughout. Use
[`pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md),
[`summarize_pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md),
and snake_case argument names such as `m_init`, `m_range`, and `q_new`.

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
