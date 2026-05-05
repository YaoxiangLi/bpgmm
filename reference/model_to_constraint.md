# Convert PGMM Paper Model Names to Constraint Codes

Convert PGMM Paper Model Names to Constraint Codes

## Usage

``` r
model_to_constraint(model)
```

## Arguments

- model:

  Character scalar naming one of the eight PGMM covariance structures:
  \`CCC\`, \`CCU\`, \`CUC\`, \`CUU\`, \`UCC\`, \`UCU\`, \`UUC\`, or
  \`UUU\`.

## Value

Integer vector of length three. \`1\` means constrained and \`0\` means
unconstrained.
