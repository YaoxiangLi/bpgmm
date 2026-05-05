# bpgmm: Bayesian inference for parsimonious Gaussian mixture models

[![CRAN
status](https://www.r-pkg.org/badges/version/bpgmm)](https://cran.r-project.org/package=bpgmm)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/bpgmm)](https://cran.r-project.org/package=bpgmm)
[![R-CMD-check](https://github.com/YaoxiangLi/bpgmm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YaoxiangLi/bpgmm/actions/workflows/R-CMD-check.yaml)
[![GitHub
release](https://img.shields.io/github/v/release/YaoxiangLi/bpgmm?display_name=tag&sort=semver)](https://github.com/YaoxiangLi/bpgmm/releases)
[![pkgdown](https://github.com/YaoxiangLi/bpgmm/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/YaoxiangLi/bpgmm/actions/workflows/pkgdown.yaml)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1007%2Fs00357--021--09391--8-blue)](https://doi.org/10.1007/s00357-021-09391-8)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

`bpgmm` implements Bayesian inference for parsimonious Gaussian mixture
models. It is designed for model-based clustering when the number of
clusters, the object partition, and the cluster covariance structure are
all inferential targets.

The package uses Markov chain Monte Carlo for posterior estimation and
reversible-jump MCMC (RJMCMC) for model selection across constrained
mixtures of factor analyzers.

## Installation

Install the released version from CRAN:

``` r

install.packages("bpgmm")
```

Install the development version from GitHub:

``` r

install.packages("remotes")
remotes::install_github("YaoxiangLi/bpgmm")
```

Load the package:

``` r

library(bpgmm)
```

## What bpgmm Does

`bpgmm` is useful when you want to:

- estimate the number of clusters from the data;
- infer posterior cluster membership probabilities;
- compare parsimonious covariance structures;
- fit mixture-of-factor-analyzers models for high-dimensional
  clustering;
- use RJMCMC to move across models with different parameter dimensions.

## Minimal Runnable Example

This example creates two small clusters, fits a short RJMCMC chain, and
summarizes the posterior samples. Real analyses should use a longer
burn-in and more posterior samples.

``` r

set.seed(2026)

X <- cbind(
  matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
  matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
)
known_labels <- rep(1:2, each = 4)

fit <- pgmm_rjmcmc(
  X = X,
  mInit = 2,
  mVec = c(1, 3),
  qnew = 1,
  burn = 1,
  niter = 3,
  constraint = model_to_constraint("UUU"),
  Mstep = 0,
  Vstep = 0
)

fit_summary <- summarize_pgmm_rjmcmc(fit, trueCluster = known_labels)
as.integer(fit_summary$nCluster["2"])
#> [1] 3

as.integer(fit_summary$nConstraint["UUU"])
#> [1] 3

fit_summary$ari
#> [1] 1
```

Here `X` is a numeric matrix with variables in rows and observations in
columns. Set `Mstep = 1` to allow RJMCMC updates for the number of
clusters and `Vstep = 1` to allow updates for the variance structure.

The main user-facing function is
[`pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md).
Important settings include:

- `mInit`: initial number of clusters.
- `mVec`: allowed cluster range, such as `c(1, 6)`.
- `qnew`: number of latent factors for a newly proposed cluster.
- `burn` and `niter`: burn-in and posterior sampling iterations.
- `constraint`: initial covariance model, usually set with
  [`model_to_constraint()`](https://yaoxiangli.github.io/bpgmm/reference/model_to_constraint.md).
- `Mstep`, `Vstep`, and `SCind`: switches for cluster-number,
  covariance-model, and split/combine RJMCMC moves.

The older names
[`pgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md),
[`summarizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md),
and the misspelled
[`summerizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md)
are deprecated compatibility wrappers. They still work, but new code
should use
[`pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md)
and
[`summarize_pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md).

The eight covariance structures in the paper are represented by model
labels such as `"CCC"`, `"CUU"`, and `"UUU"`. The package also accepts
the legacy three-number constraint encoding:

``` r

model_to_constraint("UUU")
constraint_to_model(c(1, 0, 0))
```

## Paper

The methodology behind this package is described in:

> Lu, X., Li, Y., & Love, T. (2021). On Bayesian Analysis of
> Parsimonious Gaussian Mixture Models. *Journal of Classification*, 38,
> 576-593. <https://doi.org/10.1007/s00357-021-09391-8>

The paper develops an RJMCMC inferential procedure for constrained
mixture-of-factor-analyzers models. The inferential goals are the
partition of observations, the number of clusters, and the covariance
structure of the clusters, each represented through posterior
distributions.

## Citation

If you use `bpgmm` in published work, please cite both the package and
the methodology paper. In R, run:

``` r

citation("bpgmm")
```

BibTeX for the paper:

``` bibtex
@article{lu2021bayesian,
  author = {Lu, Xiang and Li, Yaoxiang and Love, Tanzy},
  title = {On Bayesian Analysis of Parsimonious Gaussian Mixture Models},
  journal = {Journal of Classification},
  year = {2021},
  volume = {38},
  pages = {576--593},
  doi = {10.1007/s00357-021-09391-8}
}
```

## License

`bpgmm` is released under the GPL-3 license.

## Website

The package website is available at
<https://yaoxiangli.github.io/bpgmm/>. It includes a reference index,
release notes, and worked examples that run the actual package
functions.
