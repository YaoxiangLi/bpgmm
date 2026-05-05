# generatePriorPsi

generate prior value for parameter Psi

## Usage

``` r
generatePriorPsi(p, m, delta, bbeta, constraint)
```

## Arguments

- p:

  the number of features

- m:

  the number of clusters

- delta:

  hyperparameters

- bbeta:

  hyperparameters

- constraint:

  the pgmm constraint, a vector of length three with binary entry. For
  example, c(1,1,1) means the fully constraint model
