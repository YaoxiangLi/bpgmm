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
#> bpgmm 1.2.6 loaded. If you use bpgmm in published work, please cite it with citation("bpgmm").

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

``` r

cluster_cols <- c("#0072B2", "#D55E00", "#009E73")
plot(
  X[1, ], X[2, ],
  col = cluster_cols[known_labels],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Simulated data",
  asp = 1
)
legend(
  "topleft",
  legend = paste("Known", sort(unique(known_labels))),
  col = cluster_cols[sort(unique(known_labels))],
  pch = 19,
  bty = "n"
)
```

![Scatter plot of simulated observations colored by known
cluster.](examples_files/figure-html/unnamed-chunk-3-1.png)

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
For research use, increase `burn` and `niter`, and set `m_step = 1` and
`v_step = 1` when you want RJMCMC model selection for the number of
clusters and covariance structure.

``` r

fit_log <- capture.output({
  fit <- pgmm_rjmcmc(
    X = X,
    m_init = 2,
    m_range = c(1, 3),
    q_new = 1,
    burn = 1,
    niter = 3,
    constraint = constraint,
    m_step = 0,
    v_step = 0,
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
#>  [1] "tau_samples"            "psi_samples"           
#>  [3] "mean_samples"           "lambda_samples"        
#>  [5] "factor_score_samples"   "allocation_samples"    
#>  [7] "constraint_samples"     "alpha1_samples"        
#>  [9] "alpha2_samples"         "beta_samples"          
#> [11] "active_cluster_samples"
length(fit$allocation_samples)
#> [1] 3
```

## Summarize posterior samples

[`summarize_pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md)
returns the modal allocation, posterior counts for the number of
clusters, posterior counts for covariance models, and optionally an
adjusted Rand index.

``` r

fit_summary <- summarize_pgmm_rjmcmc(fit, true_cluster = known_labels)

fit_summary$allocation
#> [1] 2 2 2 2 1 1 1 1
fit_summary$n_clusters
#> 
#> 2 
#> 3
fit_summary$n_constraints
#> 
#> UUU 
#>   3
fit_summary$ari
#> [1] 1
```

``` r

old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(
  X[1, ], X[2, ],
  col = cluster_cols[known_labels],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Known labels",
  asp = 1
)
plot(
  X[1, ], X[2, ],
  col = cluster_cols[fit_summary$allocation],
  pch = 19,
  xlab = "Variable 1",
  ylab = "Variable 2",
  main = "Posterior allocation",
  asp = 1
)
```

![Two-panel plot comparing known labels and posterior
allocation.](examples_files/figure-html/unnamed-chunk-8-1.png)

``` r

par(old_par)
```

The posterior samples also give simple model-selection summaries. In
this tiny chain the settings fix the number of clusters and covariance
model, but the same code becomes more informative when `m_step = 1` and
`v_step = 1`.

``` r

old_par <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
barplot(
  fit_summary$n_clusters,
  col = "#56B4E9",
  border = NA,
  xlab = "Number of clusters",
  ylab = "Posterior sample count",
  main = "Cluster count"
)
barplot(
  fit_summary$n_constraints,
  col = "#E69F00",
  border = NA,
  xlab = "PGMM model",
  ylab = "Posterior sample count",
  main = "Covariance model"
)
```

![Posterior bar plots for cluster count and covariance
model.](examples_files/figure-html/unnamed-chunk-9-1.png)

``` r

par(old_par)
```

## Scale the example to real data

For a real analysis, keep the same workflow but use larger sampler
settings. The exact values depend on data size and convergence
diagnostics.

``` r

fit <- pgmm_rjmcmc(
  X = your_data_matrix,
  m_init = 2,
  m_range = c(1, 10),
  q_new = 4,
  burn = 5000,
  niter = 15000,
  constraint = model_to_constraint("UUU"),
  m_step = 1,
  v_step = 1,
  split_combine = 1,
  verbose = FALSE
)

summarize_pgmm_rjmcmc(fit, true_cluster = known_labels)
```
