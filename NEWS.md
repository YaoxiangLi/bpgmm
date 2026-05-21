# bpgmm 1.1.3

* Fixed `Mstep = 1` by loading the cluster-number proposal helpers as package code rather than leaving them in an ignored nested `R/` directory.
* Fixed split/combine RJMCMC moves by passing `X` explicitly to split-allocation helpers and adding the missing latent-factor update helper.
* Fixed the C++ allocation update to add `log(tao)` to component log densities instead of adding raw mixture weights.
* Added C++ input validation for multivariate normal densities, log-ratio calculations, allocation dimensions, and mixture weights.
* Expanded unit coverage for all eight PGMM covariance constraints, C++ helpers, cluster-number RJMCMC, covariance RJMCMC, and split/combine moves.

# bpgmm 1.1.2

* Modernized README with badges, installation guidance, paper citation, and model-constraint helper examples.
* Added `summarizePgmmRJMCMC()` as the correctly spelled summary function and retained `summerizePgmmRJMCMC()` for backward compatibility.
* Added tidyverse-style `pgmm_rjmcmc()` and `summarize_pgmm_rjmcmc()` as the preferred public API.
* Deprecated `pgmmRJMCMC()`, `summarizePgmmRJMCMC()`, and `summerizePgmmRJMCMC()`; these names remain available as compatibility wrappers.
* Added helpers `model_to_constraint()` and `constraint_to_model()` for translating between paper model labels and legacy constraint vectors.
* Improved package startup citation guidance for users publishing results from `bpgmm`.
* Added unit tests for the public API, covariance-constraint mapping, summary helpers, and native C++ wrappers.
* Fixed zero-iteration handling in `pgmm_rjmcmc()` and added validation for sampler inputs and summary result objects.
* Added a `verbose` argument to suppress per-iteration progress output in examples, tests, and scripted workflows.
* Updated C++ build settings and validation for current R/Rcpp best practices.
