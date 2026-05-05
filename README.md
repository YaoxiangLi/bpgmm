# bpgmm

[![CRAN status](https://www.r-pkg.org/badges/version/bpgmm)](https://cran.r-project.org/package=bpgmm)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/bpgmm)](https://cran.r-project.org/package=bpgmm)

`bpgmm` implements Bayesian inference for parsimonious Gaussian mixture
models. It is designed for model-based clustering when the number of clusters,
the object partition, and the cluster covariance structure are all inferential
targets.

The package uses Markov chain Monte Carlo for posterior estimation and
reversible-jump MCMC (RJMCMC) for model selection across constrained mixtures of
factor analyzers.

## Installation

Install the released version from CRAN:

```r
install.packages("bpgmm")
```

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("YaoxiangLi/bpgmm")
```

Load the package:

```r
library(bpgmm)
```

## What bpgmm Does

`bpgmm` is useful when you want to:

- estimate the number of clusters from the data;
- infer posterior cluster membership probabilities;
- compare parsimonious covariance structures;
- fit mixture-of-factor-analyzers models for high-dimensional clustering;
- use RJMCMC to move across models with different parameter dimensions.

The main user-facing function is `pgmmRJMCMC()`.

```r
fit <- pgmmRJMCMC(
  X = X,
  mInit = 2,
  mVec = c(1, 6),
  qnew = 2,
  burn = 100,
  niter = 1000,
  Mstep = 1,
  Vstep = 1
)
```

Here `X` is a numeric matrix with variables in rows and observations in columns.
Set `Mstep = 1` to allow RJMCMC updates for the number of clusters and
`Vstep = 1` to allow updates for the variance structure.

## Paper

The methodology behind this package is described in:

> Lu, X., Li, Y., & Love, T. (2021). On Bayesian Analysis of Parsimonious
> Gaussian Mixture Models. *Journal of Classification*, 38, 576-593.
> https://doi.org/10.1007/s00357-021-09391-8

The paper develops an RJMCMC inferential procedure for constrained
mixture-of-factor-analyzers models. The inferential goals are the partition of
observations, the number of clusters, and the covariance structure of the
clusters, each represented through posterior distributions.

## Citation

If you use `bpgmm` in published work, please cite both the package and the
methodology paper. In R, run:

```r
citation("bpgmm")
```

BibTeX for the paper:

```bibtex
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
