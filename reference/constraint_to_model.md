# Convert PGMM Constraint Codes to Paper Model Names

The paper represents the eight PGMM covariance structures with
three-letter model names. Each letter is either \`C\` for constrained or
\`U\` for unconstrained. The letters indicate whether the loading matrix
is shared across clusters, whether the noise covariance is shared across
clusters, and whether the noise covariance is isotropic within clusters.

## Usage

``` r
constraint_to_model(constraint)
```

## Arguments

- constraint:

  Integer or numeric vector of length three with entries \`0\` or \`1\`.
  \`1\` maps to \`C\`; \`0\` maps to \`U\`.

## Value

A character scalar, one of \`CCC\`, \`CCU\`, \`CUC\`, \`CUU\`, \`UCC\`,
\`UCU\`, \`UUC\`, or \`UUU\`.
