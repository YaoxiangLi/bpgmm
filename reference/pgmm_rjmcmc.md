# Bayesian Model-Based Clustering with Parsimonious Gaussian Mixture Models

Carries out model-based clustering using parsimonious Gaussian mixture
models. MCMC is used for parameter estimation and RJMCMC is used for
model selection.

## Usage

``` r
pgmm_rjmcmc(
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
  verbose = TRUE
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

## Value

A list of posterior samples with snake_case fields: \`tau_samples\`,
\`psi_samples\`, \`mean_samples\`, \`lambda_samples\`,
\`factor_score_samples\`, \`allocation_samples\`,
\`constraint_samples\`, \`alpha1_samples\`, \`alpha2_samples\`,
\`beta_samples\`, and \`active_cluster_samples\`.

## Details

The \`constraint\` argument follows the three-letter PGMM model notation
used in Lu, Li, and Love (2021). The first entry indicates whether
loading matrices are shared across clusters, the second whether noise
covariance matrices are shared across clusters, and the third whether
the noise covariance is isotropic within each cluster. Use
\[model_to_constraint()\] to convert model names such as \`CCC\`,
\`CCU\`, \`CUC\`, \`CUU\`, \`UCC\`, \`UCU\`, \`UUC\`, and \`UUU\` into
the numeric vector used internally.
