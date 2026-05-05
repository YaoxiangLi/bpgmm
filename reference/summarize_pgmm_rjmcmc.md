# Summarize RJMCMC Samples from a Bayesian PGMM Fit

Summarizes posterior samples from \[pgmm_rjmcmc()\] into the modal
allocation, posterior counts for the number of clusters, posterior
counts for the eight PGMM covariance-constraint models, and optionally
the adjusted Rand index against a known reference partition.

## Usage

``` r
summarize_pgmm_rjmcmc(pgmmResList, trueCluster = NULL)

summarizePgmmRJMCMC(pgmmResList, trueCluster = NULL)

summerizePgmmRJMCMC(pgmmResList, trueCluster = NULL)
```

## Arguments

- pgmmResList:

  Result list from \[pgmm_rjmcmc()\].

- trueCluster:

  Optional true or reference cluster allocation.

## Value

A list with \`Zalloc\`, \`nCluster\`, \`nConstraint\`, and optionally
\`ari\`.
