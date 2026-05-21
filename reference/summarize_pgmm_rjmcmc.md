# Summarize RJMCMC Samples from a Bayesian PGMM Fit

Summarizes posterior samples from \[pgmm_rjmcmc()\] into the modal
allocation, posterior counts for the number of clusters, posterior
counts for the eight PGMM covariance-constraint models, and optionally
the adjusted Rand index against a known reference partition.

## Usage

``` r
summarize_pgmm_rjmcmc(fit, true_cluster = NULL)
```

## Arguments

- fit:

  Result list from \[pgmm_rjmcmc()\].

- true_cluster:

  Optional true or reference cluster allocation.

## Value

A list with \`allocation\`, \`n_clusters\`, \`n_constraints\`, and
optionally \`ari\`.
