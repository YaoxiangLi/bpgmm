# Changelog

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
- Updated C++ build settings and validation for current R/Rcpp best
  practices.
