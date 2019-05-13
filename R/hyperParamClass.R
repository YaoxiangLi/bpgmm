#' @import methods
#' @name Hparam
#' @title Hparam-class
#' @aliases Hparam
#' @description Definiton of hyper parameter sets
#'
#' @export
setClass(
  "Hparam",
  slots = c(
    alpha1 = "numeric",
    alpha2 = "numeric",
    delta  = "numeric",
    ggamma = "numeric",
    bbeta  = "numeric"
  ),
  prototype = list(
    alpha1 = numeric(),
    alpha2 = numeric(),
    delta  = numeric(),
    ggamma = numeric(),
    bbeta  = numeric()
  )
)


#' @export
setValidity("Hparam", function(object) {
  if (object@alpha1 < 0 |
      object@alpha2 < 0 |
      object@delta  < 0 |
      object@ggamma < 0 |
      object@bbeta  < 0) {
    "Hyperparameter should be non-negative!"
  }
})
