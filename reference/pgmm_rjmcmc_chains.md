# Run Multiple Independent Bayesian PGMM Chains

Runs independent \`pgmm_rjmcmc()\` chains, optionally in parallel. This
is the safest way to use multiple CPU cores because each MCMC iteration
depends on the previous state, while independent chains can be evaluated
concurrently.

## Usage

``` r
pgmm_rjmcmc_chains(
  X,
  m_init,
  m_range,
  q_new,
  delta = 2,
  ggamma = 2,
  burn = 20,
  niter = 1000,
  constraint = c(0, 0, 0),
  d_vec = c(1, 1, 1),
  s_vec = c(1, 1, 1),
  m_step = 0,
  v_step = 0,
  split_combine = 0,
  verbose = FALSE,
  chains = 2,
  cores = min(chains, available_cores()),
  seed = NULL
)
```

## Arguments

- X:

  the observation matrix with variables in rows and observations in
  columns.

- m_init:

  the number of initial clusters.

- m_range:

  the allowed range for the number of clusters.

- q_new:

  the number of latent factors for a new cluster.

- delta:

  scalar hyperparameter for the noise covariance prior

- ggamma:

  scalar hyperparameter used in covariance-structure proposals

- burn:

  the number of burn-in iterations

- niter:

  the number of posterior sampling iterations

- constraint:

  initial PGMM covariance constraint. Use a three-letter model label
  such as \`"CCC"\` or \`"UUU"\`, or a numeric vector of length three
  with binary entries. For example, \`c(1, 1, 1)\` is \`CCC\`, the fully
  constrained model, and \`c(0, 0, 0)\` is \`UUU\`, the fully
  unconstrained model.

- d_vec:

  a vector of hyperparameters with length three, shape parameters for
  alpha1, alpha2 and bbeta respectively

- s_vec:

  a vector of hyperparameters with length three, rate parameters for
  alpha1, alpha2 and bbeta respectively

- m_step:

  indicator for RJMCMC model selection on the number of clusters.

- v_step:

  indicator for RJMCMC model selection on covariance structures.

- split_combine:

  indicator for using split/combine moves in the cluster-number RJMCMC
  step.

- verbose:

  logical; if \`TRUE\`, print iteration progress.

- chains:

  positive integer giving the number of independent chains.

- cores:

  positive integer giving the number of worker processes to use. Values
  greater than \`chains\` are reduced to \`chains\`.

- seed:

  optional integer seed used to generate deterministic per-chain seeds.

## Value

A list with one fitted \`pgmm_rjmcmc()\` result per chain. The result
has class \`bpgmm_rjmcmc_chains\` and stores the per-chain seeds in the
\`chain_seeds\` attribute.
