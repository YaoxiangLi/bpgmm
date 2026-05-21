# Changelog

## bpgmm 1.1.5

- Fixed the allocation prior contribution in RJMCMC acceptance
  calculations to use the log product of allocated mixture weights,
  matching the paper’s joint posterior.

## bpgmm 1.1.4

- Standardized internal R helper names to snake_case while preserving
  the exported compatibility wrappers.
- Added snake_case wrappers around generated Rcpp entry points and
  routed package internals through those wrappers.
- Cleaned the source layout by renaming R files to lowercase hyphenated
  names and removing the ignored duplicate `R/choosem/` tree.
- Kept result-list names and legacy public arguments stable to avoid
  breaking existing user code.

## bpgmm 1.1.3

- Fixed `Mstep = 1` by loading the cluster-number proposal helpers as
  package code rather than leaving them in an ignored nested `R/`
  directory.
- Fixed split/combine RJMCMC moves by passing `X` explicitly to
  split-allocation helpers and adding the missing latent-factor update
  helper.
- Fixed the C++ allocation update to add `log(tao)` to component log
  densities instead of adding raw mixture weights.
- Added C++ input validation for multivariate normal densities,
  log-ratio calculations, allocation dimensions, and mixture weights.
- Expanded unit coverage for all eight PGMM covariance constraints, C++
  helpers, cluster-number RJMCMC, covariance RJMCMC, and split/combine
  moves.

## bpgmm 1.1.2

- Modernized README with badges, installation guidance, paper citation,
  and model-constraint helper examples.
- Added
  [`summarizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md)
  as the correctly spelled summary function and retained
  [`summerizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md)
  for backward compatibility.
- Added tidyverse-style
  [`pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md)
  and
  [`summarize_pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md)
  as the preferred public API.
- Deprecated
  [`pgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md),
  [`summarizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md),
  and
  [`summerizePgmmRJMCMC()`](https://yaoxiangli.github.io/bpgmm/reference/summarize_pgmm_rjmcmc.md);
  these names remain available as compatibility wrappers.
- Added helpers
  [`model_to_constraint()`](https://yaoxiangli.github.io/bpgmm/reference/model_to_constraint.md)
  and
  [`constraint_to_model()`](https://yaoxiangli.github.io/bpgmm/reference/constraint_to_model.md)
  for translating between paper model labels and legacy constraint
  vectors.
- Improved package startup citation guidance for users publishing
  results from `bpgmm`.
- Added unit tests for the public API, covariance-constraint mapping,
  summary helpers, and native C++ wrappers.
- Fixed zero-iteration handling in
  [`pgmm_rjmcmc()`](https://yaoxiangli.github.io/bpgmm/reference/pgmm_rjmcmc.md)
  and added validation for sampler inputs and summary result objects.
- Added a `verbose` argument to suppress per-iteration progress output
  in examples, tests, and scripted workflows.
- Updated C++ build settings and validation for current R/Rcpp best
  practices.
