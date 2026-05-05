# bpgmm 1.1.2

* Modernized README with badges, installation guidance, paper citation, and model-constraint helper examples.
* Added `summarizePgmmRJMCMC()` as the correctly spelled summary function and retained `summerizePgmmRJMCMC()` for backward compatibility.
* Added helpers `model_to_constraint()` and `constraint_to_model()` for translating between paper model labels and legacy constraint vectors.
* Improved package startup citation guidance for users publishing results from `bpgmm`.
* Added unit tests for the public API, covariance-constraint mapping, summary helpers, and native C++ wrappers.
* Updated C++ build settings and validation for current R/Rcpp best practices.

