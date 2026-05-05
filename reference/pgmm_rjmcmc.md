# Bayesian Model-Based Clustering with Parsimonious Gaussian Mixture Models

Carries out model-based clustering using parsimonious Gaussian mixture
models. MCMC is used for parameter estimation and RJMCMC is used for
model selection.

## Usage

``` r
pgmm_rjmcmc(
  X,
  mInit,
  mVec,
  qnew,
  delta = 2,
  ggamma = 2,
  burn = 20,
  niter = 1000,
  constraint = c(0, 0, 0),
  dVec = c(1, 1, 1),
  sVec = c(1, 1, 1),
  Mstep = 0,
  Vstep = 0,
  SCind = 0
)

pgmmRJMCMC(
  X,
  mInit,
  mVec,
  qnew,
  delta = 2,
  ggamma = 2,
  burn = 20,
  niter = 1000,
  constraint = c(0, 0, 0),
  dVec = c(1, 1, 1),
  sVec = c(1, 1, 1),
  Mstep = 0,
  Vstep = 0,
  SCind = 0
)
```

## Arguments

- X:

  the observation matrix with size p \* m

- mInit:

  the number of initial clusters

- mVec:

  the range of the number of clusters

- qnew:

  the number of latent factors for a new cluster

- delta:

  scalar hyperparameter for the noise covariance prior

- ggamma:

  scalar hyperparameter used in covariance-structure proposals

- burn:

  the number of burn-in iterations

- niter:

  the number of posterior sampling iterations

- constraint:

  initial PGMM covariance constraint, a numeric vector of length three
  with binary entries. For example, \`c(1, 1, 1)\` is \`CCC\`, the fully
  constrained model, and \`c(0, 0, 0)\` is \`UUU\`, the fully
  unconstrained model.

- dVec:

  a vector of hyperparameters with length three, shape parameters for
  alpha1, alpha2 and bbeta respectively

- sVec:

  a vector of hyperparameters with length three, rate parameters for
  alpha1, alpha2 and bbeta respectively

- Mstep:

  indicator for RJMCMC model selection on the number of clusters

- Vstep:

  indicator for RJMCMC model selection on covariance structures

- SCind:

  indicator for using split/combine moves in the cluster-number RJMCMC
  step

## Details

The \`constraint\` argument follows the three-letter PGMM model notation
used in Lu, Li, and Love (2021). The first entry indicates whether
loading matrices are shared across clusters, the second whether noise
covariance matrices are shared across clusters, and the third whether
the noise covariance is isotropic within each cluster. Use
\[model_to_constraint()\] to convert model names such as \`CCC\`,
\`CCU\`, \`CUC\`, \`CUU\`, \`UCC\`, \`UCU\`, \`UUC\`, and \`UUU\` into
the numeric vector used internally.
